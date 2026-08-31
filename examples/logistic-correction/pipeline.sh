#!/usr/bin/env bash
# Logistic abundance correction vs raw counts.
#
# Default: same layout as examples/realistic (NCBI genomes + public Kraken2/Kaiju).
# After exec, compare.py writes raw vs logistic L1 and R² against true taxIDs from the
# same annotations (one run; correction is Annotation → abundance export).
#
# Reuse an existing realistic run:
#   SAMOVAR_REALISTIC_DIR=/path/to/examples/realistic/run bash pipeline.sh
# Quick toy annotators instead of public indexes:
#   SAMOVAR_TOY=1 bash pipeline.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-${SCRIPT_DIR}/run}"
realistic_dir="${SAMOVAR_REALISTIC_DIR:-${SCRIPT_DIR}/../realistic/run}"

chmod +x "${SCRIPT_DIR}/logistic_corrector.py" "${SCRIPT_DIR}/compare.py"
samovar tools import -n logistic_correction \
  --exec-path "${SCRIPT_DIR}/logistic_corrector.py" \
  --type export \
  --pytest \
  --flags "--min-efficiency 0.05" || true

reuse=""
if [[ -d "${realistic_dir}/initial_annotations" && -d "${realistic_dir}/regenerated_annotations" ]]; then
  reuse="$realistic_dir"
fi
if [[ -d "${output_dir}/initial_annotations" && -d "${output_dir}/regenerated_annotations" ]]; then
  reuse="$output_dir"
fi

if [[ -n "$reuse" && "${SAMOVAR_FORCE_RUN:-0}" != "1" ]]; then
  echo "Comparing existing run: $reuse"
  python "${SCRIPT_DIR}/compare.py" --run "$reuse" -o "${SCRIPT_DIR}/figures"
  echo "Done: ${SCRIPT_DIR}/figures"
  exit 0
fi

mkdir -p "$output_dir"

if [[ "${SAMOVAR_TOY:-0}" == "1" ]]; then
  toy_db="$(cd "${SCRIPT_DIR}/../toy" && pwd)/run/.database"
  samovar_ensure_toy_annotators "$toy_db"
  export SAMOVAR_ALLOW_TEST_GENOMES=1
  export SAMOVAR_REUSE_GENOMES=1
  samovar generate \
      --genome_dir "$SAMOVAR/data/test_genomes/meta" \
      --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
      --output_dir "$output_dir" \
      --n_samples 2 \
      --total_reads "${SAMOVAR_N_READS:-400}" \
      --host_fraction 0.2
  samovar prepare \
      --output_dir "$output_dir" \
      --test-genomes \
      --N_reads "${SAMOVAR_N_READS:-400}" \
      --export logistic \
      --kraken2-test "kraken2 toy" \
      --kaiju-test "kaiju toy"
else
  mkdir -p "$output_dir/.genomes"
  ORG_GROUPS=(Archaea Bacteria Viridiplantae Alveolata Fungi Metazoa Viruses)
  existing_genomes="$(find "${output_dir}/.genomes" -maxdepth 1 \( -name '*.fa' -o -name '*.fa.gz' -o -name '*.fna' -o -name '*.fna.gz' -o -name '*.fasta' -o -name '*.fasta.gz' \) 2>/dev/null | wc -l)"
  genome_src="${output_dir}/.genomes"
  if [[ "$existing_genomes" -lt 21 && -d "${realistic_dir}/.genomes" ]]; then
    existing_real="$(find "${realistic_dir}/.genomes" -maxdepth 1 \( -name '*.fa' -o -name '*.fa.gz' -o -name '*.fna' -o -name '*.fna.gz' -o -name '*.fasta' -o -name '*.fasta.gz' \) 2>/dev/null | wc -l)"
    if [[ "$existing_real" -ge 21 ]]; then
      echo "Reusing genomes from ${realistic_dir}/.genomes"
      genome_src="${realistic_dir}/.genomes"
    fi
  fi
  if [[ "$genome_src" == "${output_dir}/.genomes" && "$existing_genomes" -lt 21 ]]; then
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

  K2_DIR="${SAMOVAR_KRAKEN2_DB_ROOT}/standard_8GB_2025oct"
  K2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20251015.tar.gz"
  samovar_ensure_database kraken2 standard_8GB "$K2_DIR" "hash.k2d" "$K2_URL"

  KAIJU_DIR="${SAMOVAR_KAIJU_DB:-/mnt/tank/scratch/partition-metagenomics/databases/kaiju/refseq_2024aug}"
  KAIJU_URL="https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_2024-08-14.tgz"
  samovar_ensure_database kaiju refseq "$KAIJU_DIR" "*.fmi" "$KAIJU_URL"

  samovar generate \
      --genome_dir "$genome_src" \
      --host_genome "${SAMOVAR}/data/test_genomes/host/9606.fna" \
      --n_samples "${SAMOVAR_N_SAMPLES:-4}" \
      --total_reads "${SAMOVAR_N_READS:-8000}" \
      --output_dir "$output_dir" \
      --cores "${SAMOVAR_GENERATE_CORES:-4}"

  samovar prepare \
      --output_dir "$output_dir" \
      --kraken2-test "kraken2 standard_8GB" \
      --kaiju-test "kaiju refseq" \
      --max-genomes "${SAMOVAR_MAX_GENOMES:-40}" \
      --cores "${SAMOVAR_CORES:-16}" \
      --export logistic
fi

samovar_run_exec "$output_dir"
python "${SCRIPT_DIR}/compare.py" --run "$output_dir" -o "${SCRIPT_DIR}/figures"
echo "Done: $output_dir"
echo "Raw vs logistic figures: ${SCRIPT_DIR}/figures"
