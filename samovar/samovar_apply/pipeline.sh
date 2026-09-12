#!/usr/bin/env bash
# Demonstrate samovar apply and samovar apply --full on a dummy pipeline.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
# shellcheck disable=SC1091
source "$ROOT/examples/common.sh"
samovar_setup_env

DEMO="${SAMOVAR_APPLY_DEMO:-$ROOT/samovar/samovar_apply/work}"
rm -rf "$DEMO"
mkdir -p "$DEMO/reads_source" "$DEMO/reads_new" "$DEMO/source_run/.log/configs"

python3 - "$DEMO" "$ROOT" <<'PY'
import sys
from pathlib import Path

import pandas as pd
import yaml

from samovar.reprofilers import run_reprofiler

demo = Path(sys.argv[1])
root = Path(sys.argv[2])


def write_fastq(dest: Path, n: int = 8) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    r1, r2 = [], []
    for i in range(n):
        taxid = 562 if i % 2 == 0 else 9606
        rec = f"@r{i}|taxid:{taxid}\nACGTACGTACGT\n+\nIIIIIIIIIIII\n"
        r1.append(rec)
        r2.append(rec)
    (dest / "sample_R1.fastq").write_text("".join(r1))
    (dest / "sample_R2.fastq").write_text("".join(r2))


write_fastq(demo / "reads_source")
write_fastq(demo / "reads_new")
source = demo / "source_run"
n = 16
tax = [9606, 9606, 562, 562] * (n // 4)
regenerated = pd.DataFrame(
    {
        "seq": [f"t{i}" for i in range(n)],
        "taxid_dummy": tax,
        "length": [12] * n,
        "true": tax,
    }
)
initial = {
    "sample.annotation": pd.DataFrame(
        {
            "seq": ["a", "b", "c", "d"],
            "taxid_dummy": [9606, 562, 9606, 562],
            "length": [12] * 4,
        }
    )
}
ground = {"dummy": pd.DataFrame({"taxid": [9606, 562], "N_1": [8, 8]})}
regen = source / "regenerated_annotations"
regen.mkdir(parents=True, exist_ok=True)
regenerated.to_csv(regen / "combined_annotation_table.csv", index=False)
gt = source / "regenerated" / ".regenerated_abundance"
gt.mkdir(parents=True, exist_ok=True)
ground["dummy"].to_csv(gt / "dummy.csv", index=False)
run_reprofiler(
    "linear",
    regenerated=regenerated,
    ground_truth=ground,
    initial=initial,
    output_dir=source / "reprofiled_annotations",
    config={"seed": 0},
)
cfg = source / ".log" / "configs"


def dump(name, data):
    path = cfg / name
    path.write_text(yaml.safe_dump(data, sort_keys=False))


dump(
    "config_init.yaml",
    {
        "r1_dir": str(source / "initial_trimmed"),
        "r2_dir": str(source / "initial_trimmed"),
        "output_dir": str(source / "initial_reports"),
        "run_config": [
            {
                "run_name": "dummy",
                "type": "dummy",
                "cmd": "dummy",
                "db_path": ".",
                "threads": 1,
            }
        ],
    },
)
dump("config_qc.yaml", {"output_dir": str(source), "qc": ""})
dump("config_scoring.yaml", {"output_dir": str(source), "scoring_tools": []})
dump(
    "config_export.yaml",
    {"output_dir": str(source), "export": "identity", "export_formats": ["abundance"]},
)
dump(
    "config_reprofiling.yaml",
    {
        "output_dir": str(source / "reprofiled_annotations"),
        "initial_dir": str(source / "initial_annotations"),
        "regenerated_path": str(regen / "combined_annotation_table.csv"),
        "ground_truth_dir": str(gt),
        "reprofiler": "linear",
        "seed": 0,
    },
)
dump(
    "config_annotation2iss.yaml",
    {
        "annotation_dir": str(source / "initial_annotations"),
        "observed_abundance_dir": str(source / "initial_abundance"),
        "abundance_dir": str(gt),
        "output_dir": str(source / "regenerated"),
        "regeneration_mode": "direct",
        "table_reads_generator": "direct",
        "seed": 0,
        "cores": 1,
    },
)
(source / ".log" / "window.env").write_text(
    "export SAMOVAR_START=setup_reads\nexport SAMOVAR_END=viz_reprofiled\n"
)
print(f"source pipeline ready under {source}")
PY

samovar apply \
  --input_dir "$DEMO/reads_new" \
  --pipeline "$DEMO/source_run" \
  --output_dir "$DEMO/apply_normal" \
  --cores "${SAMOVAR_CORES:-1}"

samovar apply --full \
  --input_dir "$DEMO/reads_new" \
  --pipeline "$DEMO/source_run" \
  --output_dir "$DEMO/apply_full" \
  --cores "${SAMOVAR_CORES:-1}"

echo "Normal apply: $DEMO/apply_normal"
echo "Full apply:   $DEMO/apply_full"
ls -la "$DEMO/apply_normal/reprofiled_annotations" "$DEMO/apply_full/reprofiled_annotations"
