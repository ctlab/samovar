#!/usr/bin/env bash
# Realistic-scale layout: NCBI genomes + public Kraken2/Kaiju already in the catalog.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-${SCRIPT_DIR}/run}"
mkdir -p "$output_dir/.genomes"

ORG_GROUPS=(Archaea Bacteria Viridiplantae Alveolata Fungi Metazoa Viruses)
if ! samovar_seed_public_genomes "$output_dir/.genomes" 21; then
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

samovar_public_index_vars
samovar_ensure_database kraken2 "$K2_NAME" "$K2_DIR" "hash.k2d" "$K2_URL"
samovar_ensure_database kaiju "$KAIJU_NAME" "$KAIJU_DIR" "*.fmi" "$KAIJU_URL"

samovar generate \
    --genome_dir "${output_dir}/.genomes" \
    --host_genome "${SAMOVAR}/data/test_genomes/host/9606.fna" \
    --n_samples "${SAMOVAR_N_SAMPLES:-4}" \
    --total_reads "${SAMOVAR_N_READS:-8000}" \
    --output_dir "$output_dir" \
    --cores "${SAMOVAR_GENERATE_CORES:-4}"

samovar prepare \
    --output_dir "$output_dir" \
    --kraken2-test "kraken2 ${K2_NAME}" \
    --kaiju-test "kaiju ${KAIJU_NAME}" \
    --max-genomes "${SAMOVAR_MAX_GENOMES:-40}" \
    --cores "${SAMOVAR_CORES:-16}"

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
