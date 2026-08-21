#!/usr/bin/env bash
# Legacy helper. Optional MetaPhlAn bowtie2 db: METAPHLAN_BOWTIE2DB or $3.
set -euo pipefail

metaphlan_samovar () {
  local out_dir="$1"
  local read_dir="$2"
  local db="${3:-${METAPHLAN_BOWTIE2DB:-${SAMOVAR_METAPHLAN_DB:-}}}"
  local threads="${METAPHLAN_THREADS:-4}"

  local R1=($(ls -d "${read_dir}"/[1-9]_*R1* 2>/dev/null || ls -d "${read_dir}"/*R1*))
  local R2=($(ls -d "${read_dir}"/*R2*))

  mkdir -p "$out_dir"
  for i in "${!R1[@]}"; do
    local concat
    concat=$(basename "${R1[i]}" | sed 's/_.*//g')
    echo -e "\n" "$concat" "\n"

    local db_args=()
    if [[ -n "$db" ]]; then
      db_args+=(--bowtie2db "$db")
    fi

    metaphlan \
      --input_type fastq \
      --nproc "$threads" \
      "${R1[i]}" \
      --bowtie2out "$out_dir/${concat}.bowtie2.bz2" \
      -o "$out_dir/$concat.out" \
      "${db_args[@]}"
  done
}

metaphlan_samovar "${1:-reports/initial/metaphlan}" "${2:-.}" "${3:-}"
