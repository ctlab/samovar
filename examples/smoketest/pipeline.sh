#!/usr/bin/env bash
# Short generate → prepare → exec: one sample, concat regeneration, few reads.
# Toy Kraken2/Kaiju indexes are reused only when the install config already
# names them. Otherwise they are built with the same thread count as prepare.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

# Default is 3/4 of the CPUs. install.sh ask sets SAMOVAR_CORES from its prompt.
if [[ -z "${SAMOVAR_CORES:-}" ]]; then
  ncpu="$(nproc 2>/dev/null || echo 1)"
  SAMOVAR_CORES=$(( ncpu * 3 / 4 ))
fi
if (( SAMOVAR_CORES < 1 )); then
  SAMOVAR_CORES=1
fi
export SAMOVAR_CORES

output_dir="$(samovar_example_outdir)"
rm -rf "$output_dir/"

# Build a missing toy index with prepare's thread count (kaiju-mkbwt -n).
samovar_build_toy_index() {
  local kind="$1"
  local dest yaml
  dest="$(samovar_toy_database_dir)"
  mkdir -p "$dest"
  yaml="$dest/config.yaml"
  cat > "$yaml" << EOF
input_dir:
  - "${SAMOVAR}/data/test_genomes/meta"
  - "${SAMOVAR}/data/test_genomes/host"
output_dir: "${dest}/database_prep"
mutation_rate: 0.02
include_percent: 70.0
EOF
  echo "Building ${kind} toy index with ${SAMOVAR_CORES} threads"
  samovar build_database --type "$kind" --config_path "$yaml" \
    --db_path "$dest/${kind}_db" --example-omit --index toy --flags "" \
    --threads "$SAMOVAR_CORES"
}

for kind in kraken2 kaiju; do
  if [[ -z "$(samovar_catalog_database_path "$kind" toy)" ]]; then
    samovar_build_toy_index "$kind"
  else
    echo "Reusing catalog database ${kind}/toy"
  fi
done

# Two genomes, no host. A host genome starts a second ISS pass, and a third
# genome only adds stitching. ISS itself does not get faster with more CPUs
# at this read count: the error model loads on one thread, and extra workers
# write empty temp FASTQs.
samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --genomes "$SAMOVAR/data/test_genomes/meta/4932.fna" \
    --genomes "$SAMOVAR/data/test_genomes/meta/562.fna" \
    --n_samples 1 \
    --total_reads 100 \
    --cores 1 \
    --output_dir "$output_dir"

# concat merges both annotators into one abundance table, so regeneration
# writes one FASTQ. direct keeps a separate sample per annotator.
samovar prepare \
    --output_dir "$output_dir" \
    --test-genomes \
    --cores "$SAMOVAR_CORES" \
    --max-genomes 2 \
    --table_reads_generator concat \
    --N 1 \
    --N_reads 100 \
    --no-multiqc \
    --no-feature-importance \
    --endpoint reprofile \
    --kraken2-test "kraken2 toy" \
    --kaiju-test "kaiju toy"

samovar_run_exec "$output_dir"
