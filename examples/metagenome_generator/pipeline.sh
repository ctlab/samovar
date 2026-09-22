#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
rm -rf "$output_dir/"

chmod +x "${SCRIPT_DIR}/constant_iss.py"
samovar tools import -n constant_iss \
  --exec-path "${SCRIPT_DIR}/constant_iss.py" \
  --type meta \
  --flags "--model hiseq"

samovar_ensure_toy_annotators "$(samovar_toy_database_dir)"

export SAMOVAR_ALLOW_TEST_GENOMES=1
export SAMOVAR_REUSE_GENOMES=1

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --output_dir "$output_dir" \
    --n_samples 2 \
    --total_reads 200 \
    --host_fraction 0.2 \
    --metagenome_generator constant_iss \
    --flags constant_iss "--n-reads 200"

samovar prepare \
    --output_dir "$output_dir" \
    --test-genomes \
    --N_reads 200 \
    --metagenome_generator constant_iss \
    --flags metagenome_generator "--n-reads 200" \
    --kraken2-test "kraken2 toy" \
    --kaiju-test "kaiju toy"

samovar_run_exec "$output_dir"
echo "Done: $output_dir"
