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
    out = []
    for col in columns:
        name = str(col)
        name = re.sub(r"_[0-9]+$", "", name)
        name = name.replace("taxid", "taxID")
        name = re.sub(r"\.\.\.[0-9]+$", "", name)
        out.append(name)
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
        df = df.drop(df.columns[drop], axis=1)
    return df


def read_annotation_dir(
    data_dir: Union[str, Path],
    sample_name_position: int = 0,
) -> pd.DataFrame:
    """Load per-sample annotation CSVs, skipping combined tables."""
    data_dir = Path(data_dir)
    frames = []
    for path in sorted(data_dir.glob("*.csv")):
        if path.name.startswith("combined_annotation_table"):
            continue
        parts = path.name.split(".")
        sample_name = "_".join(parts[: sample_name_position + 1])
        try:
            tmp = pd.read_csv(path)
        except (pd.errors.EmptyDataError, pd.errors.ParserError, OSError):
            continue
        if tmp.empty:
            continue
        tmp.columns = _normalize_colnames(tmp.columns)
        tmp = _coalesce_duplicate_columns(tmp)
        tmp["sample"] = sample_name
        frames.append(tmp)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def annotator_columns(df: pd.DataFrame) -> list:
    cols = []
    for col in df.columns:
        name = str(col)
        if "confidence" in name.lower():
            continue
        if name.startswith("taxID_") or name.startswith("N_"):
            if name.lower() == "read_type":
                continue
            cols.append(col)
    return cols


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
