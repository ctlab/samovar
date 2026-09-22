#!/usr/bin/env bash
# Same community as examples/realistic, plus an assembly annotator (taxID and MAG_ID).
# MegaHIT/GTDB slots stay identity so the example does not need those tools.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="$(samovar_example_outdir)"
host="${SAMOVAR}/data/test_genomes/host/9606.fna"
phage_acc=(GCF_000819615.1 GCF_000840245.1 GCF_000836945.1 GCF_000844825.1)

gen_args=(
  --reindex 0
  --host_genome "$host"
  --n_samples "${SAMOVAR_N_SAMPLES:-4}"
  --total_reads "${SAMOVAR_N_READS:-8000}"
  --output_dir "$output_dir"
  --cores "${SAMOVAR_GENERATE_CORES:-4}"
)
prep_args=(
  --output_dir "$output_dir"
  --assembly-test "assembly ."
  --assembler identity
  --gene-caller identity
  --binner identity
  --binner-qc identity
  --binner-combine identity
  --mag-taxonomy identity
  --aligner identity
  --mag-quantifier identity
  --max-genomes "${SAMOVAR_MAX_GENOMES:-40}"
  --cores "${SAMOVAR_CORES:-16}"
)

if samovar_light_public_indexes; then
  SAMOVAR_PHASE=indexes bash "${SAMOVAR}/examples/phage/pipeline.sh"
  samovar generate --accessions "${phage_acc[@]}" --host_fraction 0.15 "${gen_args[@]}"
  samovar prepare \
    --kraken2-test "kraken2 phage_test" \
    --kaiju-test "kaiju phage_test" \
    "${prep_args[@]}"
else
  samovar_ensure_public_indexes
  samovar generate --accessions "${phage_acc[@]}" GCF_000005845.2 "${gen_args[@]}"
  samovar prepare \
    --kraken2-test "kraken2 ${K2_NAME}" \
    --kaiju-test "kaiju ${KAIJU_NAME}" \
    "${prep_args[@]}"
fi

samovar_run_exec "$output_dir"
samovar multiqc --output_dir "$output_dir" -- --export --interactive
samovar_harvest_example "$output_dir" "$SCRIPT_DIR"

python - <<PY
from pathlib import Path
import pandas as pd
root = Path("$output_dir")
init = next((root / "initial_annotations").glob("*.annotation.csv"))
df = pd.read_csv(init, nrows=0)
cols = list(df.columns)
assert any(str(c).startswith("feat_") and "MAG_ID" in str(c) for c in cols), cols
assert any("taxID_assembly" in str(c) or "taxid_assembly" in str(c).lower() for c in cols), cols
print("assembly MAG_ID feature ok", init)
PY
echo "Done: $output_dir"
