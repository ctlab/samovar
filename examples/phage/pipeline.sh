#!/usr/bin/env bash
# Two phage communities on named phage_test Kaiju and Kraken2 indexes.
# Kaiju includes GCF_000867865.1; Kraken2 includes GCF_000844825.1.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

part1="${SAMOVAR_OUTDIR_1:-${SAMOVAR}/examples_outdir/phage/database1}"
part2="${SAMOVAR_OUTDIR_2:-${SAMOVAR}/examples_outdir/phage/database2}"
db="${SAMOVAR}/examples_outdir/phage/database"
host="${SAMOVAR}/data/test_genomes/host/9606.fna"
n_samples="${N_SAMPLES:-2}"
total_reads="${TOTAL_READS:-2000}"

mkdir -p "$db"
if [[ ! -e "$db/kraken2_db/hash.k2d" ]] || ! find -L "$db/kaiju_db" -name '*.fmi' 2>/dev/null | grep -q .; then
  samovar genome download --output-dir "$db/kaiju_src" \
    GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000867865.1
  samovar genome download --output-dir "$db/kraken2_src" \
    GCF_000840245.1 GCF_000836945.1 GCF_000844825.1
  cat > "$db/kaiju.yaml" << EOF
input_dir:
  - ${db}/kaiju_src
output_dir: ${db}/kaiju_prep
mutation_rate: 0.0
include_percent: 100.0
EOF
  cat > "$db/kraken2.yaml" << EOF
input_dir:
  - ${db}/kraken2_src
output_dir: ${db}/kraken2_prep
mutation_rate: 0.0
include_percent: 100.0
EOF
  samovar build --type kaiju \
    --config_path "$db/kaiju.yaml" \
    --db_path "$db/kaiju_db" \
    --index phage_test --flags ""
  samovar build --type kraken2 \
    --config_path "$db/kraken2.yaml" \
    --db_path "$db/kraken2_db" \
    --index phage_test --flags ""
fi
samovar import -n phage_test --type database --tool kaiju --exec-path "$db/kaiju_db"
samovar import -n phage_test --type database --tool kraken2 --exec-path "$db/kraken2_db"

samovar generate \
  --accessions GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1 \
  --reindex 1 \
  --host_genome "$host" \
  --host_fraction 0.15 \
  --output_dir "$part1" \
  --n_samples "$n_samples" \
  --total_reads "$total_reads" \
  --cores 1
samovar prepare \
  --output_dir "$part1" \
  --kraken2-test "kraken2 phage_test" \
  --kaiju-test "kaiju phage_test"
samovar exec --output_dir "$part1"

samovar generate \
  --accessions GCF_000819615.1 GCA_035764635.1 GCF_000836945.1 GCF_000844825.1 GCF_000867865.1 \
  --reindex 0 \
  --host_genome "$host" \
  --host_fraction 0.15 \
  --output_dir "$part2" \
  --n_samples "$n_samples" \
  --total_reads "$total_reads" \
  --cores 1
samovar prepare \
  --output_dir "$part2" \
  --kraken2-test "kraken2 phage_test" \
  --kaiju-test "kaiju phage_test"
samovar exec --output_dir "$part2"

samovar multiqc --output_dir "$part1" -- --export --interactive || true
samovar_harvest_example "$part1" "$SCRIPT_DIR" || true
