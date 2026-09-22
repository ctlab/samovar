#!/usr/bin/env bash
# Logistic abundance correction vs raw counts.
#
# Default: same layout as examples/realistic (NCBI genomes + public Kraken2/Kaiju).
# After exec, compare.py writes raw vs logistic L1 and R² against true taxIDs from the
# same annotations (one run; correction is Annotation → abundance export).
#
# Reuse an existing realistic run:
#   SAMOVAR_REALISTIC_DIR=/path/to/examples_outdir/realistic bash pipeline.sh
# Quick toy annotators instead of public indexes:
#   SAMOVAR_TOY=1 bash pipeline.sh
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
realistic_dir="${SAMOVAR_REALISTIC_DIR:-${SAMOVAR}/examples_outdir/realistic}"

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
  toy_db="$(samovar_toy_database_dir)"
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
  host="${SAMOVAR}/data/test_genomes/host/9606.fna"
  phage_acc=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)
  if samovar_light_public_indexes; then
    SAMOVAR_PHASE=indexes bash "${SAMOVAR}/examples/phage/pipeline.sh"
    samovar generate \
      --accessions "${phage_acc[@]}" \
      --reindex 0 \
      --host_genome "$host" \
      --host_fraction 0.15 \
      --n_samples "${SAMOVAR_N_SAMPLES:-4}" \
      --total_reads "${SAMOVAR_N_READS:-8000}" \
      --output_dir "$output_dir" \
      --cores "${SAMOVAR_GENERATE_CORES:-4}"
    samovar prepare \
      --output_dir "$output_dir" \
      --kraken2-test "kraken2 phage_test" \
      --kaiju-test "kaiju phage_test" \
      --max-genomes "${SAMOVAR_MAX_GENOMES:-40}" \
      --cores "${SAMOVAR_CORES:-16}" \
      --export logistic
  else
    samovar_ensure_public_indexes
    samovar generate \
      --accessions "${phage_acc[@]}" GCF_000005845.2 \
      --reindex 0 \
      --host_genome "$host" \
      --n_samples "${SAMOVAR_N_SAMPLES:-4}" \
      --total_reads "${SAMOVAR_N_READS:-8000}" \
      --output_dir "$output_dir" \
      --cores "${SAMOVAR_GENERATE_CORES:-4}"
    samovar prepare \
      --output_dir "$output_dir" \
      --kraken2-test "kraken2 ${K2_NAME}" \
      --kaiju-test "kaiju ${KAIJU_NAME}" \
      --max-genomes "${SAMOVAR_MAX_GENOMES:-40}" \
      --cores "${SAMOVAR_CORES:-16}" \
      --export logistic
  fi
fi

samovar_run_exec "$output_dir"
python "${SCRIPT_DIR}/compare.py" --run "$output_dir" -o "${SCRIPT_DIR}/figures"
samovar_harvest_example "$output_dir" "$SCRIPT_DIR" || true
echo "Done: $output_dir"
echo "Raw vs logistic figures: ${SCRIPT_DIR}/figures"
