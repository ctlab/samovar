#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-samovar_out}"
mkdir -p "$output_dir/.database"

# Public Kraken2 indexes (override any URL via K2_*_URL). See:
# https://benlangmead.github.io/aws-indexes/k2
K2_MINUSB_URL="${K2_MINUSB_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_minusb_20250402.tar.gz}"
K2_STANDARD_8_URL="${K2_STANDARD_8_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20250402.tar.gz}"
K2_STANDARD_16_URL="${K2_STANDARD_16_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_standard_16gb_20250402.tar.gz}"
K2_PLUSPF_8_URL="${K2_PLUSPF_8_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_pluspf_08gb_20250402.tar.gz}"
K2_PLUSPFP_8_URL="${K2_PLUSPFP_8_URL:-https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_08gb_20250402.tar.gz}"

declare -A K2_URLS=(
  [minusb]="$K2_MINUSB_URL"
  [standard_8]="$K2_STANDARD_8_URL"
  [standard_16]="$K2_STANDARD_16_URL"
  [pluspf_8]="$K2_PLUSPF_8_URL"
  [pluspfp_8]="$K2_PLUSPFP_8_URL"
)

preprocess_args=()
for name in minusb standard_8 standard_16 pluspf_8 pluspfp_8; do
  dest="$output_dir/.database/kraken2_${name}"
  samovar_fetch_archive "${K2_URLS[$name]}" "$dest" "hash.k2d"
  preprocess_args+=(--kraken2-"${name}" "kraken2 ${dest}")
done

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --output_dir "$output_dir"

samovar preprocess \
    --output_dir "$output_dir" \
    "${preprocess_args[@]}"

samovar exec --output_dir "$output_dir"
