#!/usr/bin/env bash
# Public Kraken2 indexes vs one realistic community (NCBI genomes).
# Names and download URLs only; local runs reuse indexes already in the catalog.
# Launch with SAMOVAR_SLURM=1 SAMOVAR_SLURM_CPUS=16 for the cluster exec step.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
mkdir -p "$output_dir/.log"
host="${SAMOVAR}/data/test_genomes/host/9606.fna"
phage_acc=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)

# Catalog names + official URLs. A preinstalled index is reused from the
# SamovaR catalog; otherwise the URL is lazy-downloaded.
# pracken (~353 GB) is used only when already installed (SAMOVAR_INCLUDE_PRACKEN=1).
# SAMOVAR_CI_LIGHT_INDEXES=1 uses phage_test (built like examples/phage).
declare -A K2_URLS=(
  [standard_8GB]="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08_GB_20251015.tar.gz"
  [virus]="https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20251015.tar.gz"
  [pracken]="https://genome-idx.s3.amazonaws.com/kraken/k2_NCBI_reference_20251007.tar.gz"
)

preprocess_args=()
if samovar_light_public_indexes; then
  SAMOVAR_PHASE=indexes bash "${SAMOVAR}/examples/phage/pipeline.sh"
  preprocess_args+=(--kraken2-phage "kraken2 phage_test")
  preprocess_args+=(--kaiju-phage "kaiju phage_test")
else
  for name in standard_8GB virus; do
    resolved="$(samovar_ensure_named_database kraken2 "$name" "hash.k2d" "${K2_URLS[$name]}")"
    # CLI suffixes must not contain extra underscores (Snakemake sample wildcards).
    case "$name" in
      standard_8GB) preprocess_args+=(--kraken2-std8gb "kraken2 ${resolved}") ;;
      virus) preprocess_args+=(--kraken2-viral "kraken2 ${resolved}") ;;
    esac
  done
  if [[ "${SAMOVAR_INCLUDE_PRACKEN:-0}" == "1" ]]; then
    if resolved="$(samovar_ensure_named_database kraken2 pracken "hash.k2d" "${K2_URLS[pracken]}" 0)"; then
      preprocess_args+=(--kraken2-pracken "kraken2 ${resolved}")
    else
      echo "REPORT: pracken is not installed and is >50 GB; not downloading."
    fi
  fi
fi

gen_acc=("${phage_acc[@]}")
host_args=()
if samovar_light_public_indexes; then
  host_args=(--host_fraction 0.15)
else
  gen_acc+=(GCF_000005845.2)
fi

samovar generate \
    --accessions "${gen_acc[@]}" \
    --reindex 0 \
    --host_genome "$host" \
    "${host_args[@]}" \
    --n_samples "${SAMOVAR_N_SAMPLES:-4}" \
    --total_reads "${SAMOVAR_N_READS:-8000}" \
    --output_dir "$output_dir" \
    --cores "${SAMOVAR_GENERATE_CORES:-4}"

samovar prepare \
    --output_dir "$output_dir" \
    "${preprocess_args[@]}" \
    --cores "${SAMOVAR_CORES:-16}"

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
samovar_harvest_example "$output_dir" "$SCRIPT_DIR"
