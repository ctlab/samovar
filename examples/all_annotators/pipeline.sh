#!/usr/bin/env bash
# Four annotators on one NCBI bacterial community.
# Databases are catalog names plus their official download URLs.
# A local run reuses an index already registered in the SamovaR catalog.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
# RefSeq bacteria present in the public Kraken2 / Kaiju / KrakenUniq indexes.
bact_acc=(GCF_000005845.2 GCF_000006945.2 GCF_000009045.1)

K2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz"
KAIJU_URL="https://kaiju-idx.s3.eu-central-1.amazonaws.com/2024/kaiju_db_refseq_2024-08-14.tgz"
KRAKEN_URL="https://ccb.jhu.edu/data/minikraken/minikraken_20171013_4GB.tgz"
KUNIQ_URL="https://genome-idx.s3.amazonaws.com/kraken/kuniq_microbialdb_minus_kdb.20230808.tgz"

k2="$(samovar_ensure_named_database kraken2 standard_8GB hash.k2d "$K2_URL")"
kaiju_name="$(samovar_ensure_named_database kaiju refseq '*.fmi' "$KAIJU_URL")"
kraken_name="$(samovar_ensure_named_database kraken minikraken_4GB database.kdb "$KRAKEN_URL")"
kuniq_name="$(samovar_ensure_named_database krakenuniq microbial database.kdb "$KUNIQ_URL")"

if ! command -v kraken >/dev/null 2>&1; then
  kraken_prefix="${SAMOVAR}/examples_outdir/databases/kraken/kraken1-env"
  if [[ ! -x "${kraken_prefix}/bin/kraken" ]]; then
    echo "kraken (v1) is not on PATH; creating ${kraken_prefix}"
    conda create -y -p "$kraken_prefix" -c conda-forge -c bioconda kraken
  fi
  export PATH="${kraken_prefix}/bin:${PATH}"
fi

samovar generate \
    --accessions "${bact_acc[@]}" \
    --reindex 0 \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --n_samples "${SAMOVAR_N_SAMPLES:-3}" \
    --total_reads "${SAMOVAR_N_READS:-2000}" \
    --output_dir "$output_dir" \
    --cores "${SAMOVAR_GENERATE_CORES:-4}"

samovar prepare \
    --output_dir "$output_dir" \
    --kraken2 "kraken2 ${k2}" \
    --kaiju "kaiju ${kaiju_name}" \
    --kraken "kraken ${kraken_name}" \
    --krakenuniq "krakenuniq ${kuniq_name}" \
    --max-genomes "${SAMOVAR_MAX_GENOMES:-40}" \
    --cores "${SAMOVAR_CORES:-16}"

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
samovar_harvest_example "$output_dir" "$SCRIPT_DIR"
echo "Done: $output_dir"
