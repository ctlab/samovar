#!/usr/bin/env bash
# Legacy helper. Database path comes from $3 or KRAKEN2_DB / SAMOVAR_KRAKEN2_DB.
set -euo pipefail

kraken2_samovar () {
  local out_dir="$1"
  local read_dir="$2"
  local DB="${3:-${KRAKEN2_DB:-${SAMOVAR_KRAKEN2_DB:-}}}"
  local threads="${KRAKEN2_THREADS:-4}"

  if [[ -z "$DB" ]]; then
    echo "Set KRAKEN2_DB or pass the database directory as the third argument." >&2
    exit 1
  fi

  local R1=($(ls -d "${read_dir}"/*R1* 2>/dev/null || true))
  local R2=($(ls -d "${read_dir}"/*R2* 2>/dev/null || true))

  mkdir -p "$out_dir"
  for i in "${!R1[@]}"; do
    local concat
    concat=$(basename "${R1[i]}" | sed 's/_.*//g')
    echo -e "\n" "$concat" "\n"

    kraken2 \
      --use-names \
      --db "$DB" \
      --threads "$threads" \
      --paired "${R1[i]}" "${R2[i]}" \
      --report "$out_dir/$concat.report" \
      --output "$out_dir/$concat.out"
  done
}

kraken2_samovar "${1:-reports/initial/k2}" "${2:-.}" "${3:-}"
