#!/usr/bin/env bash
# Annotator choice and comparison: Kraken2 vs Kaiju plus an imported scoring tool.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-${SCRIPT_DIR}/run}"
rm -rf "$output_dir/"

chmod +x "${SCRIPT_DIR}/count_annotations.py"
samovar tools import -n counts \
  --exec-path "${SCRIPT_DIR}/count_annotations.py" \
  --type scoring \
  --inputs '*annotations' \
  --flags "--min-files 0"

toy_db="$(cd "${SCRIPT_DIR}/../toy" && pwd)/run/.database"
samovar_ensure_toy_annotators "$toy_db"

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --n_samples "${SAMOVAR_N_SAMPLES:-3}" \
    --total_reads "${SAMOVAR_N_READS:-400}" \
    --output_dir "$output_dir"

samovar prepare \
    --output_dir "$output_dir" \
    --test-genomes \
    --kraken2-test "kraken2 toy" \
    --kaiju-test "kaiju toy"

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
