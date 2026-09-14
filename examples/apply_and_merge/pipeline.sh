#!/usr/bin/env bash
# generate 4 samples → copy/rm into train (3) + holdout (1) → exec train →
# apply --full holdout → merge initial & regenerated → exec each merge.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=../common.sh
source "${SCRIPT_DIR}/../common.sh"
cd "$SAMOVAR"
samovar_setup_env

output_dir="${SAMOVAR_OUTDIR:-${SCRIPT_DIR}/run}"
rm -rf "$output_dir/"
mkdir -p "$output_dir/.database"
samovar_ensure_toy_annotators "$output_dir/.database"

generated="$output_dir/generated"
train="$output_dir/train"
holdout="$output_dir/holdout"
applied="$output_dir/applied"
merged_initial="$output_dir/merged_initial"
merged_regenerated="$output_dir/merged_regenerated"
n_samples="${SAMOVAR_N_SAMPLES:-4}"

samovar generate \
    --genome_dir "$SAMOVAR/data/test_genomes/meta" \
    --host_genome "$SAMOVAR/data/test_genomes/host/9606.fna" \
    --n_samples "$n_samples" \
    --total_reads "${SAMOVAR_N_READS:-400}" \
    --output_dir "$generated"

rm -rf "$train" "$holdout"
cp -a "$generated" "$train"
cp -a "$generated" "$holdout"

python3 - "$train/initial" "$holdout/initial" <<'PY'
import sys
from pathlib import Path

from samovar.seqio import find_fastq_mate, list_fastq_samples

train_dir = Path(sys.argv[1])
hold_dir = Path(sys.argv[2])
samples = list_fastq_samples(train_dir)
if len(samples) < 2:
    raise SystemExit(f"need at least 2 samples to split, got {samples!r} in {train_dir}")
hold_name = samples[-1]
train_names = samples[:-1]


def unlink_sample(folder: Path, name: str) -> None:
    for mate in ("R1", "R2"):
        path = find_fastq_mate(folder, name, mate)
        if path is not None and path.exists():
            path.unlink()


for name in train_names:
    unlink_sample(hold_dir, name)
unlink_sample(train_dir, hold_name)
left_train = list_fastq_samples(train_dir)
left_hold = list_fastq_samples(hold_dir)
print(f"train samples ({len(left_train)}): {left_train}")
print(f"holdout samples ({len(left_hold)}): {left_hold}")
if len(left_train) != len(samples) - 1 or left_hold != [hold_name]:
    raise SystemExit("split check failed")
PY

echo "directory check: train=$(python3 -c "from samovar.seqio import list_fastq_samples; print(len(list_fastq_samples('$train/initial')))") holdout=$(python3 -c "from samovar.seqio import list_fastq_samples; print(len(list_fastq_samples('$holdout/initial')))")"

samovar prepare \
    --output_dir "$train" \
    --test-genomes \
    --kraken2-test "kraken2 toy" \
    --kaiju-test "kaiju toy"

samovar_run_exec "$train"

samovar apply --full \
    --input_dir "$holdout/initial" \
    --pipeline "$train" \
    --output_dir "$applied" \
    --cores "${SAMOVAR_CORES:-1}"

samovar merge --mode initial \
    --output_dir "$merged_initial" \
    "$train" "$applied"
samovar_run_exec "$merged_initial"

samovar merge --mode regenerated \
    --output_dir "$merged_regenerated" \
    "$train" "$applied"
samovar_run_exec "$merged_regenerated"

echo "train:              $train"
echo "applied:            $applied"
echo "merged initial:     $merged_initial"
echo "merged regenerated: $merged_regenerated"
