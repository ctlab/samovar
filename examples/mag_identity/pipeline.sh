#!/usr/bin/env bash
# Toy community, identity assembly, stable MAG ids from tetranucleotide homology.
# SAMOVAR_MAG_MODE=map-initial maps regenerated reads onto the initial MAGs instead.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
rm -rf "$output_dir/"

samovar_ensure_toy_annotators "$(samovar_toy_database_dir)"

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --n_samples "${SAMOVAR_N_SAMPLES:-1}" \
    --total_reads "${SAMOVAR_N_READS:-200}" \
    --output_dir "$output_dir"

samovar prepare \
    --output_dir "$output_dir" \
    --test-genomes \
    --kraken2-test "kraken2 toy" \
    --kaiju-test "kaiju toy" \
    --assembly-test "assembly ." \
    --assembler identity \
    --gene-caller identity \
    --binner identity \
    --binner-qc identity \
    --binner-combine identity \
    --mag-taxonomy identity \
    --aligner identity \
    --mag-quantifier identity \
    --mag-mode "${SAMOVAR_MAG_MODE:-homology}" \
    --cores "${SAMOVAR_CORES:-2}"

export SAMOVAR_MAG_MODE="${SAMOVAR_MAG_MODE:-homology}"
samovar_run_exec "$output_dir"

python - <<PY
from pathlib import Path
import pandas as pd
root = Path("$output_dir")
mode = "${SAMOVAR_MAG_MODE:-homology}".replace("_", "-")
column = "MAG_ID"
init = next((root / "initial_annotations").glob("*.csv"))
df = pd.read_csv(init, nrows=0)
cols = [str(c) for c in df.columns]
assert any(c.startswith("feat_") and column in c for c in cols), cols
print("MAG feature", column, "ok", init)
PY
echo "Done: $output_dir"
