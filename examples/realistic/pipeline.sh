#!/usr/bin/env bash
# NCBI genomes + public Kraken2/Kaiju.
# SAMOVAR_CI_LIGHT_INDEXES=1 uses the phage_test community (examples/phage).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
host="${SAMOVAR}/data/test_genomes/host/9606.fna"
# Assemblies that sit in both phage_test indexes.
phage_acc=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)

if samovar_light_public_indexes; then
  # phage_test is imported by examples/phage/pipeline.sh
  samovar generate \
    --accessions "${phage_acc[@]}" \
    --reindex 0 \
    --host_genome "$host" \
    --host_fraction 0.15 \
    --n_samples "${SAMOVAR_N_SAMPLES:-4}" \
    --total_reads "${SAMOVAR_N_READS:-8000}" \
    --output_dir "$output_dir" \
    --cores "${SAMOVAR_GENERATE_CORES:-4}"
  samovar prepare \
    --output_dir "$output_dir" \
    --kraken2-test "kraken2 phage_test" \
    --kaiju-test "kaiju phage_test" \
    --max-genomes "${SAMOVAR_MAX_GENOMES:-40}" \
    --cores "${SAMOVAR_CORES:-16}"
else
  samovar_ensure_public_indexes
  samovar generate \
    --accessions "${phage_acc[@]}" GCF_000005845.2 \
    --reindex 0 \
    --host_genome "$host" \
    --n_samples "${SAMOVAR_N_SAMPLES:-4}" \
    --total_reads "${SAMOVAR_N_READS:-8000}" \
    --output_dir "$output_dir" \
    --cores "${SAMOVAR_GENERATE_CORES:-4}"
  samovar prepare \
    --output_dir "$output_dir" \
    --kraken2-test "kraken2 ${K2_NAME}" \
    --kaiju-test "kaiju ${KAIJU_NAME}" \
    --max-genomes "${SAMOVAR_MAX_GENOMES:-40}" \
    --cores "${SAMOVAR_CORES:-16}"
fi

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
samovar_harvest_example "$output_dir" "$SCRIPT_DIR"
