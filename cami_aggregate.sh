#!/usr/bin/env bash
# Collect the SAMOVAR call (column taxid_SAMOVAR) into CAMI submission files.
#
# Profiling (community percentages):
#   $OUT/reprofiled_annotations/*_reprofiled.csv
#   and, when prepare used --to-cami, the same counts already exported as
#   $OUT/exports/reprofiled/cami/SAMOVAR.profile
#   Other *.profile files in that directory are the single tools, not SAMOVAR.
#
# Taxonomy binning (one taxid per read):
#   the same reprofiled CSVs, column seq + taxid_SAMOVAR.
#
#   OUT=/mnt/san03/gentech/a.taldaev/cami/samovar/resume bash cami_aggregate.sh
#
# Writes only SAMOVAR under $OUT/cami_submit/:
#   profiling/samovar.profile
#   binning/<sample>.binning
set -euo pipefail

OUT="${OUT:-/mnt/san03/gentech/a.taldaev/cami/samovar/resume}"
DEST="${DEST:-$OUT/cami_submit}"
TAXONOMY="${TAXONOMY:-gtdb}"
TAXDUMP="${TAXDUMP:-${GTDBTK_DB:-/nvme/a.taldaev/GTDB_r232}}"

PY="${SAMOVAR_PYTHON:-}"
if [ -z "$PY" ] && [ -x /nvme/a.taldaev/miniconda3/envs/samovar/bin/python3 ]; then
    PY=/nvme/a.taldaev/miniconda3/envs/samovar/bin/python3
fi
if [ -z "$PY" ]; then
    PY="$(command -v python3)"
fi

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PYTHONPATH="${ROOT}/src${PYTHONPATH:+:${PYTHONPATH}}"
export OUT DEST TAXONOMY TAXDUMP

"$PY" - <<'PY'
import os
from collections import Counter
from pathlib import Path

import pandas as pd

from samovar.cami_profile import (
    cami_sample_id,
    counts_to_cami_flat,
    counts_to_cami_profile,
    join_cami_sections,
    try_taxonomy,
)

out = Path(os.environ["OUT"])
dest = Path(os.environ["DEST"])
taxonomy = os.environ.get("TAXONOMY") or "gtdb"
taxdump = Path(os.environ.get("TAXDUMP") or "")
src = out / "reprofiled_annotations"
if not src.is_dir():
    raise SystemExit(f"No SAMOVAR tables in {src}")

skip = {"combined_annotation_table.csv", "taxon_efficiency.json"}
frames = []
for path in sorted(src.glob("*.csv")):
    if path.name in skip or path.name.startswith("."):
        continue
    frame = pd.read_csv(path)
    if "taxid_SAMOVAR" not in frame.columns:
        continue
    if "seq" not in frame.columns:
        raise SystemExit(f"{path.name} has taxid_SAMOVAR but no seq column")
    if "sample" not in frame.columns:
        stem = path.name
        for suffix in ("_reprofiled.csv", ".csv"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
                break
        frame = frame.copy()
        frame["sample"] = stem
    frames.append(frame[["seq", "sample", "taxid_SAMOVAR"]])

if not frames:
    raise SystemExit(f"No taxid_SAMOVAR tables in {src}")

table = pd.concat(frames, ignore_index=True)
table["taxid_SAMOVAR"] = pd.to_numeric(table["taxid_SAMOVAR"], errors="coerce").fillna(0).astype(int)
table["seq"] = table["seq"].astype(str)
table["sample"] = table["sample"].map(cami_sample_id)

tax = None
if taxdump:
    for candidate in (taxdump, taxdump / "taxdump", taxdump / "taxonomy"):
        if (candidate / "nodes.dmp").is_file():
            tax = try_taxonomy(candidate, taxonomy=taxonomy)
            if tax is not None:
                print(f"[cami] taxonomy {taxonomy} from {candidate}")
                break
if tax is None:
    print("[cami] no nodes.dmp; profiling TAXPATH is the taxid itself")

profile_parts = []
bin_dir = dest / "binning"
bin_dir.mkdir(parents=True, exist_ok=True)
for sample, group in table.groupby("sample", sort=True):
    counts = Counter(int(t) for t in group["taxid_SAMOVAR"])
    if tax is not None:
        profile_parts.append(counts_to_cami_profile(counts, tax, sample))
    else:
        profile_parts.append(counts_to_cami_flat(counts, sample, taxonomy=taxonomy))
    lines = [
        "@Version:0.9.1",
        f"@SampleID:{sample}",
        f"@TaxonomyID:{taxonomy}",
        "@@SEQUENCEID\tTAXID",
    ]
    for seq, taxid in zip(group["seq"], group["taxid_SAMOVAR"]):
        lines.append(f"{seq}\t{int(taxid)}")
    path = bin_dir / f"{sample}.binning"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[cami] binning {path} ({len(group)} reads)")

prof_dir = dest / "profiling"
prof_dir.mkdir(parents=True, exist_ok=True)
profile = prof_dir / "samovar.profile"
profile.write_text(join_cami_sections(profile_parts), encoding="utf-8")
print(f"[cami] profiling {profile} ({table['sample'].nunique()} samples)")
PY
