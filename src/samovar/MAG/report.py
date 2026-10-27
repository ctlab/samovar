"""MAG mapping losses and MAG_ID shares for the MultiQC report.

Loss fractions come from ``*.mapping.json``. MAG_ID percentages are read
counts. Tables larger than a few megabytes are subsampled.
"""
from __future__ import annotations

import csv
import json
import random
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

LOSS_KEYS = ("unmapped", "multi_hit", "unresolved", "unclassified")
MAX_ROWS = 8000
FULL_SCAN_BYTES = 8_000_000


def subsample_lines(
    path: Path,
    max_rows: int = MAX_ROWS,
    full_scan_bytes: int = FULL_SCAN_BYTES,
    seed: int = 0,
) -> Tuple[str, List[str]]:
    """Return the header and at most ``max_rows`` data lines.

    Small files are read through. Larger files are sampled by byte offset so
    the cost stays bounded by the subsample, not the read count.
    """
    size = path.stat().st_size
    with path.open("rb") as handle:
        header = handle.readline().decode("utf-8", errors="replace")
        if size <= full_scan_bytes:
            body = [line.decode("utf-8", errors="replace") for line in handle]
            if len(body) <= max_rows:
                return header, body
            step = len(body) / max_rows
            return header, [body[int(i * step)] for i in range(max_rows)]
        start = handle.tell()
        span = max(size - start, 1)
        rng = random.Random(seed)
        picked: List[str] = []
        for _ in range(max_rows):
            handle.seek(start + rng.randrange(span))
            handle.readline()
            line = handle.readline()
            if line:
                picked.append(line.decode("utf-8", errors="replace"))
        return header, picked


def _delimiter(header: str) -> str:
    return "\t" if "\t" in header else ","


def _cells(line: str, delimiter: str) -> List[str]:
    return next(csv.reader([line], delimiter=delimiter))


def _mag_column_indexes(header: Sequence[str]) -> List[int]:
    indexes = []
    for i, name in enumerate(header):
        low = name.strip().lower()
        if low == "mag_id" or low.endswith("_mag_id"):
            indexes.append(i)
    return indexes


def _column_label(name: str) -> str:
    low = name.strip()
    if low.lower() == "mag_id":
        return ""
    if low.lower().endswith("_mag_id"):
        return low[: -len("_MAG_ID")]
    return low


def mag_id_shares(
    path: Path,
    *,
    max_rows: int = MAX_ROWS,
    full_scan_bytes: int = FULL_SCAN_BYTES,
) -> List[Tuple[str, Dict[str, float]]]:
    """Percent of reads per MAG_ID. One entry per MAG column."""
    header, lines = subsample_lines(path, max_rows=max_rows, full_scan_bytes=full_scan_bytes)
    if not header or not lines:
        return []
    delimiter = _delimiter(header)
    names = _cells(header, delimiter)
    indexes = _mag_column_indexes(names)
    if not indexes:
        return []
    counts = {index: {} for index in indexes}
    seen = {index: 0 for index in indexes}
    for line in lines:
        if not line.strip():
            continue
        cells = _cells(line, delimiter)
        for index in indexes:
            if index >= len(cells):
                continue
            mag = cells[index].strip() or "(unassigned)"
            bucket = counts[index]
            bucket[mag] = bucket.get(mag, 0) + 1
            seen[index] += 1
    shares = []
    for index in indexes:
        total = seen[index]
        if not total:
            continue
        label = _column_label(names[index])
        shares.append((label, {mag: round(100.0 * n / total, 2) for mag, n in counts[index].items()}))
    return shares


def loss_fractions(payload: dict) -> Dict[str, float]:
    total = float(payload.get("total") or 0)
    out = {}
    for key in LOSS_KEYS:
        fraction = payload.get(f"{key}_fraction")
        if fraction is None and total:
            fraction = float(payload.get(key) or 0) / total
        out[key] = round(100.0 * float(fraction or 0), 2)
    return out


def _known_samples(folder: Path) -> List[str]:
    if not folder.is_dir():
        return []
    suffix = ".annotation.csv"
    return [path.name[: -len(suffix)] for path in folder.glob(f"*{suffix}") if path.is_file()]


def _label(path: Path, samples: Sequence[str]) -> str:
    stem = path.name
    if stem.endswith(".mapping.json"):
        stem = stem[: -len(".mapping.json")]
    for suffix in (".assembly.out", ".annotation.csv"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    for sample in sorted(samples, key=len, reverse=True):
        if stem == sample:
            return sample
        prefix = sample + "_"
        if stem.startswith(prefix):
            rest = stem[len(prefix) :]
            return f"{sample} {rest}".strip()
    return stem


def _write_bars(series: Dict[str, Dict[str, float]], path: Path, **kwargs) -> None:
    if not series:
        return
    from samovar.stage_report import write_bargraph_mqc

    write_bargraph_mqc(series, path, ymin=0, ymax=100, xlab="Percent of reads", **kwargs)


def write_mag_mapping_plots(root: Path) -> List[Path]:
    """Write per-phase loss and MAG_ID bar graphs into the plot folders."""
    root = Path(root)
    phases = (
        ("initial", "initial_reports", "initial_annotations", "initial_annotations_plots", "Raw"),
        ("regenerated", "regenerated_reports", "regenerated_annotations", "regenerated_annotations_plots", "Regenerated"),
        ("reprofiled", "", "reprofiled_annotations", "reprofiled_annotations_plots", "Reprofiled"),
    )
    written: List[Path] = []
    for phase, reports_name, annotations_name, plots_name, title in phases:
        reports = root / reports_name if reports_name else None
        annotations = root / annotations_name
        samples = _known_samples(annotations)
        losses: Dict[str, Dict[str, float]] = {}
        if reports is not None and reports.is_dir():
            for path in sorted(reports.glob("*.mapping.json")):
                try:
                    payload = json.loads(path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    continue
                if isinstance(payload, dict) and payload.get("total"):
                    losses[_label(path, samples)] = loss_fractions(payload)
        compositions: Dict[str, Dict[str, float]] = {}
        sources: Iterable[Path] = []
        if reports is not None and reports.is_dir():
            sources = sorted(reports.glob("*.assembly.out"))
        if not sources and annotations.is_dir():
            sources = sorted(annotations.glob("*.annotation.csv"))
        for path in sources:
            for column, shares in mag_id_shares(path):
                label = _label(path, samples)
                if column:
                    label = f"{label} {column}".strip()
                compositions[label] = shares
                if phase == "reprofiled" and label not in losses:
                    losses[label] = {"unassigned": shares.get("(unassigned)", 0.0)}
        if not losses and not compositions:
            continue
        plots = root / plots_name
        plots.mkdir(parents=True, exist_ok=True)
        if losses:
            dest = plots / f"mag_{phase}_losses_mqc.json"
            _write_bars(
                losses,
                dest,
                section_name=f"{title} — MAG mapping losses",
                description=(
                    "Percent of reads. unmapped, multi-hit, unresolved, and unclassified "
                    "come from the mapping sidecar. Reprofiled bars without a sidecar "
                    "count reads whose MAG_ID is empty."
                ),
            )
            written.append(dest)
        if compositions:
            dest = plots / f"mag_{phase}_ids_mqc.json"
            _write_bars(
                compositions,
                dest,
                section_name=f"{title} — MAG_ID read shares",
                description=(
                    "Percent of reads assigned to each MAG_ID, including unassigned. "
                    f"Large tables use at most {MAX_ROWS} reads."
                ),
            )
            written.append(dest)
    return written
