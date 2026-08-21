#!/usr/bin/env bash
# Legacy helper. Database path comes from $3 or KAIJU_DB / SAMOVAR_KAIJU_DB.
set -euo pipefail

kaiju_samovar () {
  local out_dir="$1"
  local read_dir="$2"
  local DB="${3:-${KAIJU_DB:-${SAMOVAR_KAIJU_DB:-}}}"
  local threads="${KAIJU_THREADS:-4}"

  if [[ -z "$DB" ]]; then
    echo "Set KAIJU_DB or pass the database directory as the third argument." >&2
    exit 1
  fi

  local fmi
  fmi="$(find -L "$DB" -name '*.fmi' -print -quit)"
  if [[ -z "$fmi" ]]; then
    echo "No .fmi file under ${DB}" >&2
    exit 1
  fi

  local R1=($(ls -d "${read_dir}"/*R1* 2>/dev/null || true))
  local R2=($(ls -d "${read_dir}"/*R2* 2>/dev/null || true))

  mkdir -p "$out_dir"
  for i in "${!R1[@]}"; do
    local concat
    concat=$(basename "${R1[i]}" | sed 's/_.*//g')
    echo -e "\n" "$concat" "\n"

    kaiju \
      -t "$DB/nodes.dmp" \
      -f "$fmi" \
      -i "${R1[i]}" \
      -j "${R2[i]}" \
      -z "$threads" \
      -o "$out_dir/$concat.raw"

    kaiju-addTaxonNames \
      -t "$DB/nodes.dmp" \
      -n "$DB/names.dmp" \
      -i "$out_dir/$concat.raw" \
      -o "$out_dir/$concat.out"
  done
}

kaiju_samovar "${1:-reports/initial/kaiju}" "${2:-.}" "${3:-}"
