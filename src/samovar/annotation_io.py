"""Python analog of R ``read_annotation_dir`` / ``annotation2samovar``.

Used by visualization and optional abundance-table regeneration so prepare/exec
do not depend on the R package.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Iterable, Optional, Union

import pandas as pd


def _normalize_colnames(columns: Iterable) -> list:
    """Canonicalize taxID spelling; keep ``taxID_<tool>_<n>`` when several DBs share a tool."""
    prepared = []
    for col in columns:
        name = str(col)
        name = name.replace("taxid", "taxID")
        name = re.sub(r"\.\.\.[0-9]+$", "", name)
        stripped = re.sub(r"_[0-9]+$", "", name)
        prepared.append((name, stripped))
    tax_stripped = [s for name, s in prepared if name.startswith("taxID_")]
    collisions = {s for s in tax_stripped if tax_stripped.count(s) > 1}
    out = []
    for name, stripped in prepared:
        if name.startswith("taxID_") and stripped in collisions:
            out.append(name)
        elif name.startswith("taxID_") or name.startswith("N_"):
            out.append(stripped if name.startswith("taxID_") else name)
        elif name.startswith("feat_"):
            out.append(name)
        else:
            out.append(stripped)
    return out


def _coalesce_duplicate_columns(df: pd.DataFrame) -> pd.DataFrame:
    seen = {}
    drop = []
    for i, col in enumerate(df.columns):
        if col not in seen:
            seen[col] = i
            continue
        first = seen[col]
        df.iloc[:, first] = df.iloc[:, first].combine_first(df.iloc[:, i])
        drop.append(i)
    if drop:
        keep = [i for i in range(df.shape[1]) if i not in set(drop)]
        df = df.iloc[:, keep]
    return df


def contract_column_names(columns) -> list:
    """taxID / feat columns that the reprofiler treats as one feature set."""
    names = []
    for col in columns:
        name = str(col)
        low = name.lower()
        if low.startswith("taxid_") or name.startswith("feat_"):
            names.append(name)
    return names


def align_contract_columns(*frames: pd.DataFrame) -> list:
    """Give every frame the same tax and feature columns.

    A tool that classified only on one side (or returned taxid 0) still occupies
    its column on the other side, filled with 0 for tax IDs and blank for features.
    """
    import warnings

    ordered = []
    seen = set()
    for frame in frames:
        if not isinstance(frame, pd.DataFrame):
            continue
        for name in contract_column_names(frame.columns):
            if name not in seen:
                seen.add(name)
                ordered.append(name)
    aligned = []
    for frame in frames:
        if not isinstance(frame, pd.DataFrame):
            aligned.append(frame)
            continue
        missing = [name for name in ordered if name not in frame.columns]
        if not missing:
            aligned.append(frame)
            continue
        work = frame.copy()
        for name in missing:
            work[name] = 0 if name.lower().startswith("taxid_") else ""
        warnings.warn(
            "Combined annotation columns differ ("
            + ", ".join(missing)
            + "). Missing columns are filled so the ML feature set matches.",
            UserWarning,
            stacklevel=2,
        )
        aligned.append(work)
    return aligned


def _read_annotation_csv(path: Path) -> Optional[pd.DataFrame]:
    try:
        tmp = pd.read_csv(path)
    except (pd.errors.EmptyDataError, pd.errors.ParserError, OSError):
        return None
    if tmp.empty:
        return None
    tmp.columns = _normalize_colnames(tmp.columns)
    return _coalesce_duplicate_columns(tmp)


def read_annotation_dir(
    data_dir: Union[str, Path],
    sample_name_position: int = 0,
) -> pd.DataFrame:
    """Load per-sample annotation CSVs.

    When those files are gone, read ``combined_annotation_table*.csv`` and keep
    its ``sample`` column.
    """
    data_dir = Path(data_dir)
    frames = []
    for path in sorted(data_dir.glob("*.csv")):
        if path.name.startswith("combined_annotation_table"):
            continue
        parts = path.name.split(".")
        sample_name = "_".join(parts[: sample_name_position + 1])
        tmp = _read_annotation_csv(path)
        if tmp is None:
            continue
        tmp["sample"] = sample_name
        frames.append(tmp)
    if not frames:
        for path in sorted(data_dir.glob("combined_annotation_table*.csv")):
            tmp = _read_annotation_csv(path)
            if tmp is None:
                continue
            if "sample" not in tmp.columns:
                tmp["sample"] = path.stem
            frames.append(tmp)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def annotator_columns(df: pd.DataFrame) -> list:
    from samovar.annotation_columns import tax_annotator_columns

    cols = tax_annotator_columns(df.columns, include_n=True)
    return [c for c in cols if str(c).lower() != "read_type"]


def annotation_to_abundance(
    data: pd.DataFrame,
    n_reads: Optional[int] = None,
) -> Dict[str, pd.DataFrame]:
    """Count taxIDs per sample for each annotator (Python ``annotation2samovar``).

    Returns ``{annotator: DataFrame(taxid, N_<sample>...)}``. When ``n_reads``
    is set, counts are rescaled so each sample sums to that many reads.
    """
    result: Dict[str, pd.DataFrame] = {}
    if data is None or data.empty:
        return result
    sample_col = "sample" if "sample" in data.columns else None
    for col in annotator_columns(data):
        work = data[[col] + ([sample_col] if sample_col else [])].copy()
        if sample_col is None:
            work["sample"] = "1"
            sample_col_use = "sample"
        else:
            sample_col_use = sample_col
        work[col] = work[col].astype(str)
        counts = (
            work.groupby([col, sample_col_use], dropna=False)
            .size()
            .reset_index(name="n")
        )
        wide = counts.pivot(index=col, columns=sample_col_use, values="n").fillna(0)
        if n_reads:
            totals = wide.sum(axis=0).replace(0, 1)
            wide = (wide / totals * int(n_reads)).round()
        wide = wide.reset_index()
        wide.columns = ["taxid"] + [f"N_{c}" for c in wide.columns[1:]]
        result[str(col)] = wide
    return result
