"""Annotation / Feature column contract.

One tool may emit many taxonomic ``taxID_`` columns and/or many ``feat_``
columns. Scoring and abundance use only tax columns. ML/GLM use both tax
votes and ``feat_`` values.
"""

from __future__ import annotations

import re
from typing import Iterable, List, Optional, Sequence, Tuple

import pandas as pd

META_COLUMNS = frozenset(
    {"seq", "sample", "true", "read_type", "length", "classified", "taxa"}
)
DROP_RAW_COLUMNS = frozenset(
    {
        "classified",
        "taxa",
        "k-mer",
        "kmer",
        "k_mer",
        "sample",
        "true",
        "read_type",
    }
)
_TAX_PREFIX = re.compile(r"^(taxid|tax)([_.].+)?$", re.I)
_FEAT_PREFIX = re.compile(r"^feat([_.].+)?$", re.I)
_CONF = re.compile(r"confidence|_conf$", re.I)


def _name(col) -> str:
    return str(col).strip()


def is_meta_column(col) -> bool:
    return _name(col).lower() in META_COLUMNS


def is_feat_column(col) -> bool:
    """True for combined ``feat_`` columns (never scored as taxonomy)."""
    name = _name(col)
    if not name or _CONF.search(name):
        return False
    return bool(_FEAT_PREFIX.match(name))


def is_tax_column(col, *, include_true: bool = False) -> bool:
    """True for taxonomic annotation columns (``taxID_`` / ``taxid_`` / ``tax_``).

    ``feat_`` never matches. ``taxa`` (Kraken name field) is not a tax column.
    """
    name = _name(col)
    low = name.lower()
    if not name or _CONF.search(name):
        return False
    if is_feat_column(name):
        return False
    if low == "true":
        return bool(include_true)
    if low in {"taxa", "taxon", "taxonomy"}:
        return False
    if low == "taxid" or low.startswith("taxid_") or low.startswith("taxid."):
        return True
    if low.startswith("tax_") :
        return True
    return False


def tax_annotator_columns(columns: Iterable, *, include_n: bool = False) -> List:
    """Columns scored as annotators: ``taxID_*`` (and optional ``N_*``)."""
    out = []
    for col in columns:
        name = _name(col)
        if is_tax_column(name):
            out.append(col)
            continue
        if include_n and name.startswith("N_"):
            out.append(col)
    return out


def feat_annotator_columns(columns: Iterable) -> List:
    return [col for col in columns if is_feat_column(col)]


def sanitize_column_id(raw: str) -> str:
    text = re.sub(r"[^A-Za-z0-9]+", "_", str(raw).strip())
    text = text.strip("_")
    return text


def _strip_known_prefix(name: str, prefixes: Sequence[str]) -> str:
    low = name.lower()
    for prefix in prefixes:
        if low == prefix.lower():
            return ""
        token = prefix.lower().rstrip("_") + "_"
        if low.startswith(token):
            return name[len(token) :]
    return name


def tax_column_id(raw: str) -> str:
    return sanitize_column_id(
        _strip_known_prefix(_name(raw), ("taxID", "taxid", "tax_id", "tax"))
    )


def feat_column_id(raw: str) -> str:
    name = _name(raw)
    if _FEAT_PREFIX.match(name):
        rest = _strip_known_prefix(name, ("feat",))
        return sanitize_column_id(rest)
    return sanitize_column_id(name)


def classify_raw_columns(columns: Iterable) -> Tuple[Optional[object], List, List]:
    """Split a raw tool table into seq / tax / feature columns."""
    seq_col = None
    tax_cols: List = []
    feat_cols: List = []
    for col in columns:
        name = _name(col)
        low = name.lower()
        if low in {"seq", "read_id", "readid", "sequenceid", "anonymous_read_id"}:
            if seq_col is None:
                seq_col = col
            continue
        if low in DROP_RAW_COLUMNS or low in {"k-mer", "kmer"}:
            continue
        if is_tax_column(name) or low in {"taxid", "tax_id"}:
            tax_cols.append(col)
            continue
        if is_feat_column(name) or low not in {"seq", "sample", "true", "read_type"}:
            feat_cols.append(col)
    return seq_col, tax_cols, feat_cols


def combined_column_name(prefix: str, tool: str, run_id: int, col_id: str, n_same: int) -> str:
    """``taxID_tool_0`` / ``feat_tool_0_AA``. ``col_id`` may be empty when n_same==1."""
    base = f"{prefix}_{tool}_{run_id}"
    if col_id and n_same > 1:
        return f"{base}_{col_id}"
    if col_id and n_same == 1:
        return f"{base}_{col_id}" if prefix == "feat" or col_id else base
    return base


def prefix_tool_columns(df: pd.DataFrame, tool: str, run_id: int) -> pd.DataFrame:
    """Rename tax/feature columns of one tool into combined-table form.

    Tax: ``taxID_{tool}_{run_id}[_id]``. Feature: ``feat_{tool}_{run_id}[_id]``.
    Feature ``id`` is omitted only when the extractor emits a single unnamed column.
    """
    if df is None or df.empty:
        return pd.DataFrame()
    work = df.copy()
    seq_col, tax_cols, feat_cols = classify_raw_columns(work.columns)
    if seq_col is not None and seq_col != "seq":
        work = work.rename(columns={seq_col: "seq"})
        seq_col = "seq"
    if "seq" not in work.columns:
        if work.index.name == "seq" or (work.index.name is None and seq_col is None):
            work = work.reset_index()
            if "index" in work.columns and "seq" not in work.columns:
                work = work.rename(columns={"index": "seq"})
        seq_col, tax_cols, feat_cols = classify_raw_columns(work.columns)

    keep = []
    rename = {}
    tax_ids = [tax_column_id(c) for c in tax_cols]
    feat_ids = [feat_column_id(c) for c in feat_cols]
    n_tax = len(tax_cols)
    n_feat = len(feat_cols)
    for col, cid in zip(tax_cols, tax_ids):
        rename[col] = combined_column_name("taxID", tool, run_id, cid, n_tax)
        keep.append(col)
    for col, cid in zip(feat_cols, feat_ids):
        if n_feat == 1 and not cid:
            rename[col] = f"feat_{tool}_{run_id}"
        else:
            rename[col] = combined_column_name("feat", tool, run_id, cid or str(len(rename)), n_feat)
        keep.append(col)
    cols = []
    if "seq" in work.columns:
        cols.append("seq")
    cols.extend(keep)
    out = work.loc[:, [c for c in cols if c in work.columns]].rename(columns=rename)
    return out


def looks_like_header_row(fields: Sequence[str]) -> bool:
    if not fields:
        return False
    first = str(fields[0]).strip().lower().lstrip("#@")
    return first in {"seq", "read_id", "readid", "sequenceid", "anonymous_read_id"}


def select_scoring_annotators(work: Optional[pd.DataFrame], names: Optional[Sequence] = None) -> List:
    """Annotator scoring uses tax columns only, never ``feat_``.

    After prefix-strip (``taxID_kaiju_0`` → ``kaiju``) names have no prefix;
    ``feat_`` columns are still dropped by prefix.
    """
    if work is None:
        return []
    cols = list(work.columns)
    prefixed = tax_annotator_columns(cols)
    meta = {c for c in cols if is_meta_column(c) or str(c).lower() in META_COLUMNS}
    feat = {c for c in cols if is_feat_column(c)}
    if prefixed:
        pool = list(prefixed)
    else:
        pool = [c for c in cols if c not in meta and c not in feat]
    if names:
        wanted = {str(n) for n in names}
        matched = [c for c in pool if str(c) in wanted]
        if matched:
            pool = matched
        else:
            pool = [c for c in pool if not is_feat_column(c)]
    return [c for c in pool if not is_feat_column(c) and str(c).lower() != "read_type"]
