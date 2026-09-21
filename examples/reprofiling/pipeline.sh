#!/usr/bin/env bash
# Linear ML reprofiler on toy Kraken2/Kaiju labels.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
rm -rf "$output_dir/"

chmod +x "${SCRIPT_DIR}/linear_classifier.py"
samovar tools import -n linear \
  --exec-path "${SCRIPT_DIR}/linear_classifier.py" \
  --type ml \
  --flags "--max-iter 500"

toy_db="$(samovar_toy_database_dir)"
samovar_ensure_toy_annotators "$toy_db"

export SAMOVAR_ALLOW_TEST_GENOMES=1
export SAMOVAR_REUSE_GENOMES=1

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --output_dir "$output_dir" \
    --n_samples 2 \
    --total_reads 200 \
    --host_fraction 0.2

samovar prepare \
    --output_dir "$output_dir" \
    --test-genomes \
    --N_reads 200 \
    --reprofiler linear \
    --flags linear "--max-iter 400" \
    --kraken2-test "kraken2 toy" \
    --kaiju-test "kaiju toy"

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
samovar_harvest_example "$output_dir" "$SCRIPT_DIR"
echo "Done: $output_dir"
ls -l "$output_dir/reprofiled_annotations/trained_model.joblib" \
  "$output_dir/reprofiled_annotations"/*_reprofiled.csv
