"""Empty / unclassified annotator output: warn, skip regen, fail if nothing left."""

from __future__ import annotations

import logging
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd

logger = logging.getLogger(__name__)

UNCLASSIFIED = frozenset(
    {"0", "unclassified", "unclassified_root", "nan", "none", "na", ""}
)

PathLike = Union[str, Path]


class EmptyAnnotatorsError(RuntimeError):
    """Every annotator produced empty or unclassified-only assignments."""


def is_unclassified_taxid(value: Any) -> bool:
    token = str(value or "").strip().lower()
    if token.endswith(".0") and token[:-2].isdigit():
        token = token[:-2]
    return token in UNCLASSIFIED or token in {"<na>", "null"}


def _n_sample_columns(frame: pd.DataFrame) -> List[str]:
    return [c for c in frame.columns if str(c).startswith("N_")]


def abundance_has_classified_taxa(frame: pd.DataFrame) -> bool:
    if frame is None or frame.empty or "taxid" not in frame.columns:
        return False
    n_cols = _n_sample_columns(frame)
    if not n_cols:
        return False
    work = frame.copy()
    work["taxid"] = work["taxid"].astype(str)
    counts = work[n_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
    keep = (~work["taxid"].map(is_unclassified_taxid)) & (counts.sum(axis=1) > 0)
    return bool(keep.any())


def annotation_column_has_classified(frame: pd.DataFrame, col: str) -> bool:
    if frame is None or frame.empty or col not in frame.columns:
        return False
    series = frame[col].dropna().astype(str)
    return any(not is_unclassified_taxid(v) for v in series)


def _report_files_for_annotator(reports_dir: Path, name: str) -> List[Path]:
    if not reports_dir.is_dir():
        return []
    key = str(name).strip()
    if not key:
        return []
    found: List[Path] = []
    for path in sorted(reports_dir.iterdir()):
        if not path.is_file():
            continue
        low = path.name.lower()
        if not low.endswith(".out"):
            continue
        if f"_{key}." in low or f"_{key}." in path.name or f".custom_{key}." in low:
            found.append(path)
            continue
        # ``1_kaiju.kaiju.out`` / ``1_kaiju-test.kaiju.out``
        stem = path.name[: -len(".out")]
        if stem.endswith(f"_{key}") or stem.endswith(f".{key}"):
            found.append(path)
    return found


def diagnose_report_files(paths: Sequence[Path]) -> str:
    """``did_not_run`` | ``unclassified`` | ``classified``."""
    if not paths:
        return "did_not_run"
    saw_row = False
    classified = False
    nonempty = False
    for path in paths:
        try:
            if not path.is_file() or path.stat().st_size == 0:
                continue
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        lines = [ln for ln in text.splitlines() if ln.strip()]
        if not lines:
            continue
        nonempty = True
        for line in lines:
            parts = line.split("\t")
            if not parts:
                continue
            flag = parts[0].strip().upper()
            if flag in {"C", "U"} and len(parts) >= 3:
                saw_row = True
                tax = parts[2].strip()
                if "(" in tax and "taxid" in tax.lower():
                    import re

                    match = re.search(r"taxid\s*([0-9]+)", tax, re.I)
                    tax = match.group(1) if match else tax
                if flag == "C" and not is_unclassified_taxid(tax):
                    classified = True
                continue
            if len(parts) >= 2:
                saw_row = True
                if not is_unclassified_taxid(parts[1]):
                    classified = True
    if classified:
        return "classified"
    if saw_row or nonempty:
        return "unclassified"
    return "did_not_run"


def diagnose_annotator(
    name: str,
    frame: Optional[pd.DataFrame] = None,
    *,
    reports_dir: Optional[PathLike] = None,
) -> str:
    """Why this annotator has no usable taxa: did_not_run vs unclassified."""
    reports = _report_files_for_annotator(Path(reports_dir), name) if reports_dir else []
    from_reports = diagnose_report_files(reports) if reports else None
    table_ok = abundance_has_classified_taxa(frame) if frame is not None else False
    if table_ok:
        return "classified"
    if from_reports in {"did_not_run", "unclassified", "classified"}:
        if from_reports == "classified":
            return "unclassified"
        return from_reports
    if frame is None or getattr(frame, "empty", True):
        return "did_not_run"
    return "unclassified"


def warning_message(name: str, diagnosis: str) -> str:
    base = (
        f"Annotator {name!r} produced empty or unclassified-only output; "
        "it will be skipped for table regeneration."
    )
    if diagnosis == "did_not_run":
        detail = (
            " The classifier did not write output (missing or empty .out): "
            "it likely failed to start, or reads were empty and the tool was skipped."
        )
    else:
        detail = (
            " The classifier ran but assigned no taxa in this community "
            "(all labels unclassified / taxid 0)."
        )
    return base + detail


def filter_classified_abundance_tables(
    tables: Dict[str, pd.DataFrame],
    *,
    reports_dir: Optional[PathLike] = None,
    fatal_if_none: bool = True,
) -> Dict[str, pd.DataFrame]:
    """Drop empty/unclassified annotators. Raise if none remain and ``fatal_if_none``."""
    kept: Dict[str, pd.DataFrame] = {}
    dropped: List[Tuple[str, str]] = []
    for name, frame in (tables or {}).items():
        if abundance_has_classified_taxa(frame):
            kept[str(name)] = frame
            continue
        diagnosis = diagnose_annotator(str(name), frame, reports_dir=reports_dir)
        dropped.append((str(name), diagnosis))
        msg = warning_message(str(name), diagnosis)
        logger.warning(msg)
        warnings.warn(msg, UserWarning, stacklevel=2)
    if kept:
        return kept
    if not (tables or {}):
        return {}
    names = ", ".join(n for n, _ in dropped) or "all annotators"
    msg = (
        f"All annotators produced empty or unclassified-only output ({names}); "
        "cannot regenerate tables."
    )
    if fatal_if_none:
        raise EmptyAnnotatorsError(msg)
    logger.error(msg)
    return {}
