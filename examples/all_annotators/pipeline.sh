#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-samovar_10bac}"
mkdir -p "$output_dir/.database" "$output_dir/.genomes"

# Get genomes from NCBI (not from a local store)
python -m samovar.genome_fetcher \
    --output-dir "$output_dir/.genomes" \
    --N 10 \
    --group "Bacteria" \
    --email "$NCBI_EMAIL" \
    --silent

# Databases: env overrides, else build from downloaded genomes + package host FASTA
cat > "$output_dir/.database/config.yaml" << EOF
input_dir:
  - "${output_dir}/.genomes"
  - "${SAMOVAR}/data/test_genomes/host"
output_dir: "${output_dir}/database_prep"
mutation_rate: 0.02
include_percent: 70.0
EOF

if [[ -n "${SAMOVAR_KRAKEN2_DB:-}" ]]; then
    DB_KRAKEN2="$SAMOVAR_KRAKEN2_DB"
else
    DB_KRAKEN2="$output_dir/.database/kraken2_db"
    samovar build_database --type kraken2 \
        --config_path "$output_dir/.database/config.yaml" \
        --db_path "$DB_KRAKEN2"
fi

if [[ -n "${SAMOVAR_KAIJU_DB:-}" ]]; then
    DB_KAIJU="$SAMOVAR_KAIJU_DB"
else
    DB_KAIJU="$output_dir/.database/kaiju_db"
    samovar build_database --type kaiju \
        --config_path "$output_dir/.database/config.yaml" \
        --db_path "$DB_KAIJU"
fi

if [[ -n "${SAMOVAR_KRAKENUNIQ_DB:-}" ]]; then
    DB_KRAKENUNIQ="$SAMOVAR_KRAKENUNIQ_DB"
else
    DB_KRAKENUNIQ="$output_dir/.database/krakenuniq_db"
    samovar build_database --type krakenunique \
        --config_path "$output_dir/.database/config.yaml" \
        --db_path "$DB_KRAKENUNIQ"
fi

# Kraken 1 has no in-package builder; download MiniKraken unless SAMOVAR_KRAKEN_DB is set.
MINIKRAKEN_URL="${MINIKRAKEN_URL:-https://ccb.jhu.edu/software/kraken/dl/minikraken_20171019_4GB.tgz}"
if [[ -n "${SAMOVAR_KRAKEN_DB:-}" ]]; then
    DB_KRAKEN="$SAMOVAR_KRAKEN_DB"
else
    DB_KRAKEN="$output_dir/.database/kraken_db"
    samovar_fetch_archive "$MINIKRAKEN_URL" "$DB_KRAKEN" "database.kdb"
fi

samovar generate \
    --genome_dir "$output_dir/.genomes" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --n_samples 3 \
    --output_dir "$output_dir"

samovar preprocess \
    --output_dir "$output_dir" \
    --kraken2 "kraken2 $DB_KRAKEN2" \
    --kaiju "kaiju $DB_KAIJU" \
    --kraken "kraken $DB_KRAKEN" \
    --krakenuniq "krakenuniq $DB_KRAKENUNIQ"

samovar exec --output_dir "$output_dir"
