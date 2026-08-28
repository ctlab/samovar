#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-${SAMOVAR}/samovar/reprofiling_linear}"
rm -rf "$output_dir/"
mkdir -p "$output_dir/.database"

chmod +x "${SCRIPT_DIR}/linear_classifier.py"
samovar tools import -n linear \
  --exec-path "${SCRIPT_DIR}/linear_classifier.py" \
  --type ml \
  --flags "--max-iter 500"

if [[ -f "${SCRIPT_DIR}/../metagenome_generator/constant_iss.py" ]]; then
  chmod +x "${SCRIPT_DIR}/../metagenome_generator/constant_iss.py"
  samovar tools import -n constant_iss \
    --exec-path "${SCRIPT_DIR}/../metagenome_generator/constant_iss.py" \
    --type meta \
    --flags "--model hiseq"
fi

cat > "$output_dir/.database/config.yaml" << EOF
input_dir:
  - "${SAMOVAR}/data/test_genomes/meta"
  - "${SAMOVAR}/data/test_genomes/host"
output_dir: "${output_dir}/database_prep"
mutation_rate: 0.02
include_percent: 70.0
EOF

if [[ -n "${SAMOVAR_KRAKEN2_DB:-}" && -d "$SAMOVAR_KRAKEN2_DB" ]]; then
  ln -sfn "$SAMOVAR_KRAKEN2_DB" "$output_dir/.database/kraken2_db"
elif [[ -d "${SAMOVAR}/samovar_toy/.database/kraken2_db" ]] && ls "${SAMOVAR}/samovar_toy/.database/kraken2_db"/*.k2d >/dev/null 2>&1; then
  ln -sfn "${SAMOVAR}/samovar_toy/.database/kraken2_db" "$output_dir/.database/kraken2_db"
elif [[ -d "${SAMOVAR}/samovar/metagenome_generator/.database/kraken2_db" ]]; then
  ln -sfn "${SAMOVAR}/samovar/metagenome_generator/.database/kraken2_db" "$output_dir/.database/kraken2_db"
elif [[ ! -e "$output_dir/.database/kraken2_db/hash.k2d" && ! -e "$output_dir/.database/kraken2_db/taxo.k2d" ]]; then
  samovar build_database --type kraken2 --config_path "$output_dir/.database/config.yaml" --db_path "$output_dir/.database/kraken2_db" --example-omit
fi
if [[ -n "${SAMOVAR_KAIJU_DB:-}" && -e "$SAMOVAR_KAIJU_DB" ]]; then
  ln -sfn "$SAMOVAR_KAIJU_DB" "$output_dir/.database/kaiju_db"
elif [[ -f "${SAMOVAR}/samovar_toy/.database/kaiju_db/kaiju_db.fmi" ]]; then
  ln -sfn "${SAMOVAR}/samovar_toy/.database/kaiju_db" "$output_dir/.database/kaiju_db"
elif [[ -d "${SAMOVAR}/samovar/metagenome_generator/.database/kaiju_db" ]]; then
  ln -sfn "${SAMOVAR}/samovar/metagenome_generator/.database/kaiju_db" "$output_dir/.database/kaiju_db"
elif ! find "$output_dir/.database/kaiju_db" -name '*.fmi' 2>/dev/null | grep -q .; then
  samovar build_database --type kaiju --config_path "$output_dir/.database/config.yaml" --db_path "$output_dir/.database/kaiju_db" --example-omit
fi

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
    --reprofiler linear \
    --flags linear "--max-iter 400" \
    --kraken2-test "kraken2 $output_dir/.database/kraken2_db" \
    --kaiju-test "kaiju $output_dir/.database/kaiju_db"

samovar exec --output_dir "$output_dir"
echo "Done: $output_dir"
ls -l "$output_dir/reprofiled_annotations/trained_model.joblib" \
  "$output_dir/reprofiled_annotations"/*_reprofiled.csv
