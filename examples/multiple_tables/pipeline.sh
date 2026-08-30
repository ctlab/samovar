#!/usr/bin/env bash
# Several table regenerators scored with sparsedossa2-cv (table choice).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-${SCRIPT_DIR}/run}"
rm -rf "$output_dir/"

toy_db="$(cd "${SCRIPT_DIR}/../toy" && pwd)/run/.database"
samovar_ensure_toy_annotators "$toy_db"

export SAMOVAR_ALLOW_TEST_GENOMES=1
export SAMOVAR_REUSE_GENOMES=1

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --output_dir "$output_dir" \
    --n_samples 6 \
    --total_reads 400 \
    --host_fraction 0.2

samovar prepare \
    --output_dir "$output_dir" \
    --test-genomes \
    --N_reads 200 \
    --N 6 \
    --table_reads_generator direct \
    --table_reads_generator bootstrap \
    --table_reads_generator glm \
    --table_reads_generator sparsedossa2-fit \
    --table_reads_generator sparsedossa2-stool \
    --table-score sparsedossa2-cv \
    --flags sparsedossa2 "--workers ${SAMOVAR_SD2_WORKERS:-2} --cv-folds 2 --lambdas 1 --maxit 5 --max-eval 8 --prec-bits 53 --verbose" \
    --kraken2-test "kraken2 toy" \
    --kaiju-test "kaiju toy"

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
if [[ -f "$output_dir/regenerated/.regenerated_abundance/table_selection.json" ]]; then
  echo "Table selection:"
  cat "$output_dir/regenerated/.regenerated_abundance/table_selection.json"
fi
echo "Done: $output_dir"
