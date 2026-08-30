#!/usr/bin/env bash
# Public Kraken2 indexes vs one realistic community (NCBI genomes).
# Indexes live under SAMOVAR_KRAKEN2_DB_ROOT (override on other machines).
# Launch with SAMOVAR_SLURM=1 SAMOVAR_SLURM_CPUS=16 for the cluster exec step.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-${SCRIPT_DIR}/run}"
mkdir -p "$output_dir/.genomes" "$output_dir/.log"

# Same genome panel as examples/realistic (3 assemblies per NCBI group).
ORG_GROUPS=(Archaea Bacteria Viridiplantae Alveolata Fungi Metazoa Viruses)
existing_genomes="$(find "${output_dir}/.genomes" -maxdepth 1 \( -name '*.fa' -o -name '*.fa.gz' -o -name '*.fna' -o -name '*.fna.gz' -o -name '*.fasta' -o -name '*.fasta.gz' \) 2>/dev/null | wc -l)"
if [[ "$existing_genomes" -ge 21 ]]; then
  echo "Found ${existing_genomes} genomes under ${output_dir}/.genomes; skipping fetch"
else
  for org_group in "${ORG_GROUPS[@]}"; do
    tmp_dir="${output_dir}/.genomes/_tmp_${org_group}"
    rm -rf "$tmp_dir"
    mkdir -p "$tmp_dir"
    echo "Fetching 3 genomes for ${org_group}..."
    python -m samovar.genome_fetcher \
      --output-dir "$tmp_dir" \
      --N 3 \
      --group "$org_group" \
      --max-genome-mb 100 \
      --email "$NCBI_EMAIL" \
      --silent
    sleep 2
    shopt -s nullglob
    for f in "${tmp_dir}"/*-processed.fasta "${tmp_dir}"/*-processed.fasta.gz \
             "${tmp_dir}"/*.fa.gz "${tmp_dir}"/*.fna.gz "${tmp_dir}"/*.fasta.gz \
             "${tmp_dir}"/*.fa "${tmp_dir}"/*.fna "${tmp_dir}"/*.fasta; do
      [[ -f "$f" ]] || continue
      base="$(basename "$f")"
      dest="${output_dir}/.genomes/${base}"
      if [[ -e "$dest" ]]; then
        dest="${output_dir}/.genomes/${org_group}_${base}"
      fi
      mv "$f" "$dest"
    done
    shopt -u nullglob
    rm -rf "$tmp_dir"
  done
fi

K2_ROOT="${SAMOVAR_KRAKEN2_DB_ROOT}"
# Catalog names already on this cluster. pracken (~353 GB) is imported if present
# but skipped by default (SAMOVAR_INCLUDE_PRACKEN=1 to annotate with it).
declare -A K2_URLS=(
  [standard_8GB]="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20251015.tar.gz"
  [virus]="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz"
  [pracken]="https://genome-idx.s3.amazonaws.com/kraken/k2_NCBI_reference_20251007.tar.gz"
)
declare -A K2_DIRS=(
  [standard_8GB]="${K2_ROOT}/standard_8GB_2025oct"
  [virus]="${K2_ROOT}/virus_2025oct"
  [pracken]="${K2_ROOT}/pracken_2025oct"
)

preprocess_args=()
for name in standard_8GB virus; do
  samovar_ensure_database kraken2 "$name" "${K2_DIRS[$name]}" "hash.k2d" "${K2_URLS[$name]}"
done
# CLI suffixes must not contain extra underscores (Snakemake sample wildcards).
preprocess_args+=(--kraken2-std8gb "kraken2 standard_8GB")
preprocess_args+=(--kraken2-viral "kraken2 virus")
if [[ "${SAMOVAR_INCLUDE_PRACKEN:-0}" == "1" ]]; then
  if [[ -f "${K2_DIRS[pracken]}/hash.k2d" ]]; then
    samovar_ensure_database kraken2 pracken "${K2_DIRS[pracken]}" "hash.k2d" "${K2_URLS[pracken]}"
    preprocess_args+=(--kraken2-pracken "kraken2 pracken")
  else
    echo "REPORT: pracken is not on disk and is >50 GB; not downloading. Set SAMOVAR_INCLUDE_PRACKEN=1 after placing it under ${K2_DIRS[pracken]}"
  fi
fi

samovar generate \
    --genome_dir "$output_dir/.genomes" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --n_samples "${SAMOVAR_N_SAMPLES:-4}" \
    --total_reads "${SAMOVAR_N_READS:-8000}" \
    --output_dir "$output_dir" \
    --cores "${SAMOVAR_GENERATE_CORES:-4}"

samovar prepare \
    --output_dir "$output_dir" \
    "${preprocess_args[@]}" \
    --cores "${SAMOVAR_CORES:-16}"

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
