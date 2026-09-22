#!/usr/bin/env bash
# Basic SamovaR run on bundled test genomes and toy Kraken2/Kaiju indexes.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
rm -rf "$output_dir/"

samovar_ensure_toy_annotators "$(samovar_toy_database_dir)"

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --output_dir "$output_dir"

samovar prepare \
    --output_dir "$output_dir" \
    --test-genomes \
    --kraken2-test "kraken2 toy" \
    --kaiju-test "kaiju toy"

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
samovar_harvest_example "$output_dir" "$SCRIPT_DIR"
