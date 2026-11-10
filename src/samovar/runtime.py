"""Elapsed time for a checkpoint and for one annotator or profiler sample.

Absolute time is the wall clock of a checkpoint, not the sum of its jobs.
Per-sample time is the Snakemake benchmark of that sample. A missing file is
not a measured zero.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional


def record_stage(output_dir: str, stage: str, started: str, finished: str) -> Path:
    begin = float(started)
    end = float(finished)
    elapsed = end - begin
    if elapsed < 0:
        raise ValueError(f"stage {stage} finished before it started")
    dest = Path(output_dir).expanduser() / "runtime" / "stages" / f"{stage}.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(
        json.dumps(
            {"stage": stage, "started": begin, "finished": end, "elapsed_s": elapsed},
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    return dest


def _elapsed_seconds(path: Path) -> Optional[float]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if len(lines) < 2:
        return None
    header = lines[0].split("\t")
    if "s" not in header:
        return None
    try:
        value = float(lines[1].split("\t")[header.index("s")])
    except (ValueError, IndexError):
        return None
    if value < 0:
        return None
    return value


def load_sample_runtimes(root: Path) -> List[Dict[str, Any]]:
    """One row per sample job. ``run_name`` is the annotator or profiler identity."""
    rows: List[Dict[str, Any]] = []
    for phase, base in (
        ("initial", root / "initial_reports" / "runtime" / "samples"),
        ("regenerated", root / "regenerated_reports" / "runtime" / "samples"),
    ):
        if not base.is_dir():
            continue
        for path in sorted(base.glob("*/*/*.tsv")):
            elapsed = _elapsed_seconds(path)
            if elapsed is None:
                continue
            rows.append(
                {
                    "phase": phase,
                    "run_name": path.parts[-3],
                    "rule": path.parts[-2],
                    "sample": path.stem,
                    "elapsed_s": elapsed,
                }
            )
    return rows


def load_stage_runtimes(root: Path) -> List[Dict[str, Any]]:
    folder = root / "runtime" / "stages"
    if not folder.is_dir():
        return []
    rows = []
    for path in sorted(folder.glob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            elapsed = float(payload["elapsed_s"])
        except (OSError, json.JSONDecodeError, KeyError, TypeError, ValueError):
            continue
        if elapsed < 0:
            continue
        rows.append({"stage": str(payload.get("stage") or path.stem), "elapsed_s": elapsed})
    return rows


def mean_per_run(rows: Iterable[Dict[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[tuple, List[float]] = {}
    for row in rows:
        grouped.setdefault((row["phase"], row["run_name"]), []).append(float(row["elapsed_s"]))
    out = []
    for (phase, run_name), values in sorted(grouped.items()):
        out.append(
            {
                "phase": phase,
                "run_name": run_name,
                "n": len(values),
                "mean_elapsed_s": sum(values) / len(values),
            }
        )
    return out


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="samovar.runtime")
    sub = parser.add_subparsers(dest="cmd", required=True)
    record = sub.add_parser("record-stage")
    record.add_argument("output_dir")
    record.add_argument("stage")
    record.add_argument("started")
    record.add_argument("finished")
    args = parser.parse_args(argv)
    if args.cmd == "record-stage":
        record_stage(args.output_dir, args.stage, args.started, args.finished)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
