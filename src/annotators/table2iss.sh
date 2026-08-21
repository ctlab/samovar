#!/usr/bin/env bash
# Legacy ISS helper. Genome dir and output dir come from env or arguments.
set -euo pipefail

path="${1:-${SAMOVAR_GENOME_DIR:-data/test_genomes}}"
out_root="${2:-${SAMOVAR_ISS_OUT:-tests_outs/benchmarking/initial}}"
n_reads="${ISS_N_READS:-}"

mkdir -p "$out_root"
cd "$out_root"

iss_generate () {
  iss generate \
    --genomes "$1" \
    --model hiseq \
    --output "${2}_cont" \
    --n_reads "$3"
}

iss_generate_all () {
    local sample="$1"
    local human_reads="$2"
    local meta_reads="$3"
    for i in "$path"/*.fna; do
        local n="$meta_reads"
        if [[ "$(basename "$i")" == *9606* ]]; then
            n="$human_reads"
        fi
        iss_generate "$i" "$(basename "$i" .fna)_${sample}" "${n:-1000}"
    done
}

for ((i=1;i<=25;i++)); do
  human_reads=$((($RANDOM % 100) * 250))
  meta_reads=$(((25000 - $human_reads)/3))

  iss_generate_all "$i" "$human_reads" "$meta_reads"

  cat *${i}*_R1* > ${i}_full_R1.fastq
  cat *${i}*_R2* > ${i}_full_R2.fastq

  rm -f *tmp* *_cont_* *abundance* *vcf

  echo -e "done " $i "\n"
done
