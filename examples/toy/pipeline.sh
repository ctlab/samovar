#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-samovar_toy}"
rm -rf "$output_dir/"
mkdir -p "$output_dir/.database"

# Build indexes from genomes shipped with the package (no cluster paths).
# TEST/TOY EXAMPLE ONLY: Kraken2 omits Escherichia; Kaiju omits Phage Phi X.
# That gap is declared as a warning during build — do not use this for real DBs.
cat > "$output_dir/.database/config.yaml" << EOF
input_dir:
  - "${SAMOVAR}/data/test_genomes/meta"
  - "${SAMOVAR}/data/test_genomes/host"
output_dir: "${output_dir}/database_prep"
mutation_rate: 0.02
include_percent: 70.0
EOF

samovar build_database --type kraken2 --config_path "$output_dir/.database/config.yaml" --db_path "$output_dir/.database/kraken2_db" --example-omit
samovar build_database --type kaiju --config_path "$output_dir/.database/config.yaml" --db_path "$output_dir/.database/kaiju_db" --example-omit

# Optional public Kaiju index instead of the toy build:
# samovar_fetch_archive \
#   "${KAIJU_INDEX_URL:-https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_fungi_2024-08-16.tgz}" \
#   "$output_dir/.database/kaiju_db" "*.fmi"

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --output_dir "$output_dir"

samovar prepare \
    --output_dir "$output_dir" \
    --test-genomes \
    --kraken2-test "kraken2 $output_dir/.database/kraken2_db" \
    --kaiju-test "kaiju $output_dir/.database/kaiju_db"

samovar exec --output_dir "$output_dir"
# exec already runs MultiQC when it is installed; this re-runs with --export.
samovar multiqc --output_dir "$output_dir" -- --export --interactive
