"""Abundance-table contract for table regenerate and ISS.

Minimal regenerate I/O is ``taxid`` + ``N_<sample>`` columns:

    abundance table  ->  regenerated abundance table

Long annotation tables (``taxID_*`` / single-tool ``taxID``) convert to that
shape by counting reads. Phyloseq-style OTU tables (taxa × named samples,
no ``N_`` prefix) are accepted and rewritten to ``N_<sample>``.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import pandas as pd

PathLike = Union[str, Path]

INITIAL_ABUNDANCE_DIR = "initial_abundance"
REGENERATED_ABUNDANCE_DIR = Path("regenerated") / ".regenerated_abundance"

_ID_ALIASES = (
    "taxid",
    "tax_id",
    "otu",
    "otu_id",
    "#otu id",
    "feature",
    "featureid",
    "feature_id",
    "taxa",
    "taxon",
    "asv",
)
_META_COLS = frozenset(
    {
        "seq",
        "length",
        "sample",
        "true",
        "read_type",
        "taxonomy",
        "rank",
        "lineage",
        "confidence",
    }
)


def n_sample_columns(frame: pd.DataFrame) -> List[str]:
    return [c for c in frame.columns if str(c).startswith("N_")]


def _id_column(frame: pd.DataFrame) -> Optional[str]:
    lower = {str(c).strip().lower(): c for c in frame.columns}
    for alias in _ID_ALIASES:
        if alias in lower:
            return lower[alias]
    return None


def looks_like_annotation(frame: pd.DataFrame) -> bool:
    if frame is None or frame.empty:
        return False
    cols = [str(c) for c in frame.columns]
    if any(c.startswith("taxID_") or c.startswith("taxid_") for c in cols):
        return True
    if "seq" in {c.lower() for c in cols} and any(
        c.lower() in {"taxid", "taxid"} or c.lower().startswith("taxid") for c in cols
    ):
        # Long per-read table, not taxa×sample counts.
        if not n_sample_columns(frame) and "seq" in {c.lower() for c in cols}:
            return True
    return False


def is_abundance_table(frame: pd.DataFrame) -> bool:
    """True if ``frame`` is taxa × samples (already counts, not per-read calls)."""
    if frame is None or getattr(frame, "empty", True):
        return False
    if looks_like_annotation(frame) and not n_sample_columns(frame):
        return False
    if n_sample_columns(frame) and _id_column(frame):
        return True
    id_col = _id_column(frame)
    if id_col is None:
        return False
    numeric = [
        c
        for c in frame.columns
        if c != id_col and str(c).lower() not in _META_COLS and _is_numeric_series(frame[c])
    ]
    return len(numeric) >= 1


def _is_numeric_series(series: pd.Series) -> bool:
    if pd.api.types.is_numeric_dtype(series):
        return True
    coerced = pd.to_numeric(series, errors="coerce")
    return bool(coerced.notna().mean() >= 0.8)


def normalize_abundance_table(frame: pd.DataFrame) -> pd.DataFrame:
    """Return ``taxid`` plus ``N_<sample>`` columns (numeric, zeros filled)."""
    if frame is None or frame.empty:
        return pd.DataFrame(columns=["taxid"])
    work = frame.copy()
    id_col = _id_column(work)
    if id_col is None:
        work = work.reset_index()
        id_col = _id_column(work) or work.columns[0]
    n_cols = n_sample_columns(work)
    if not n_cols:
        n_cols = [
            c
            for c in work.columns
            if c != id_col and str(c).lower() not in _META_COLS and _is_numeric_series(work[c])
        ]
        renamed = {}
        for col in n_cols:
            sample = str(col)
            if sample.lower().startswith("n_"):
                sample = sample[2:]
            sample = re.sub(r"[^A-Za-z0-9._-]+", "_", sample).strip("_") or "1"
            renamed[col] = f"N_{sample}"
        work = work.rename(columns=renamed)
        n_cols = list(renamed.values())
    out = pd.DataFrame({"taxid": work[id_col].astype(str)})
    for col in n_cols:
        out[col] = pd.to_numeric(work[col], errors="coerce").fillna(0.0)
    if n_cols:
        out = out[out[n_cols].sum(axis=1) > 0]
    return out.reset_index(drop=True)


def abundance_to_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    table = normalize_abundance_table(frame)
    n_cols = n_sample_columns(table)
    if not n_cols:
        return pd.DataFrame()
    mat = table.set_index("taxid")[n_cols]
    mat.columns = [str(c)[2:] if str(c).startswith("N_") else str(c) for c in mat.columns]
    return mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)


def taxid_annotation_columns(frame: pd.DataFrame) -> List[str]:
    """Per-read taxID columns on a long annotation table."""
    cols: List[str] = []
    for col in frame.columns:
        name = str(col)
        low = name.lower()
        if name in {"seq", "length", "sample", "true", "read_type"}:
            continue
        if "confidence" in low or name.endswith("_conf"):
            continue
        if low in _ID_ALIASES and n_sample_columns(frame):
            continue
        if name.startswith("taxID") or name.startswith("taxid") or low == "taxid":
            cols.append(col)
    return cols


def annotator_name(col: str) -> str:
    name = re.sub(r"^taxID_?", "", str(col), flags=re.I)
    name = re.sub(r"^taxid_?", "", name, flags=re.I)
    name = re.sub(r"_[0-9]+$", "", name)
    if not name or name.lower() in {"taxid", "taxid"}:
        return "annotator"
    return name


def annotation_to_abundance_tables(frame: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Count per-read taxIDs into ``{annotator: taxid × N_sample}`` tables."""
    if frame is None or frame.empty:
        return {}
    work = frame
    if "sample" not in work.columns:
        work = work.copy()
        work["sample"] = "1"
    result: Dict[str, pd.DataFrame] = {}
    for col in taxid_annotation_columns(work):
        piece = work[[col, "sample"]].copy()
        piece[col] = piece[col].astype(str)
        mat = piece.groupby([col, "sample"], dropna=False).size().unstack(fill_value=0)
        mat.index = mat.index.astype(str)
        out = mat.reset_index().rename(columns={col: "taxid"})
        out.columns = ["taxid"] + [f"N_{c}" for c in out.columns[1:]]
        n_cols = n_sample_columns(out)
        if n_cols:
            out = out[out[n_cols].sum(axis=1) > 0]
        result[annotator_name(col)] = out.reset_index(drop=True)
    return result


def input_to_abundance_tables(
    data: Any,
    default_name: str = "table",
) -> Dict[str, pd.DataFrame]:
    """Accept Annotation, long annotation DataFrame, abundance table, or dict."""
    stored = getattr(data, "abundance_tables", None)
    if stored:
        return {str(k): normalize_abundance_table(v) for k, v in stored.items() if v is not None}

    frame = getattr(data, "DataFrame", data)
    if isinstance(data, dict) and not isinstance(data, pd.DataFrame):
        out: Dict[str, pd.DataFrame] = {}
        for key, value in data.items():
            if isinstance(value, pd.DataFrame):
                if is_abundance_table(value):
                    out[str(key)] = normalize_abundance_table(value)
                elif looks_like_annotation(value):
                    out.update(annotation_to_abundance_tables(value))
        return out
    if frame is None:
        return {}
    if not isinstance(frame, pd.DataFrame):
        return {}
    if frame.empty:
        return {}
    if is_abundance_table(frame):
        return {default_name: normalize_abundance_table(frame)}
    return annotation_to_abundance_tables(frame)


def dir_looks_like_annotation(path: PathLike) -> bool:
    folder = Path(path)
    if not folder.is_dir():
        return False
    for csv in sorted(folder.glob("*.csv")):
        if csv.name.startswith("combined_annotation_table") or csv.name.startswith("."):
            continue
        try:
            peek = pd.read_csv(csv, nrows=5)
        except Exception:
            continue
        if looks_like_annotation(peek) and not is_abundance_table(peek):
            return True
    return False


def load_table_input(
    source: PathLike,
    data: Optional[Any] = None,
    default_name: str = "table",
) -> Any:
    """Load annotation dir, abundance CSV, or abundance dir; or return ``data``."""
    from samovar.annotation_io import read_annotation_dir
    from samovar.parse_annotators import Annotation

    if data is not None:
        if isinstance(data, Annotation):
            return data
        tables = input_to_abundance_tables(data, default_name=default_name)
        if tables and (not isinstance(data, pd.DataFrame) or is_abundance_table(data)):
            return Annotation.from_abundance_tables(tables)
        if isinstance(data, pd.DataFrame):
            return Annotation.from_long_table(data)
        return Annotation.from_abundance_tables(tables)

    path = Path(source)
    if path.is_file():
        frame = pd.read_csv(path)
        if is_abundance_table(frame):
            return Annotation.from_abundance_tables({path.stem: frame})
        return Annotation.from_long_table(frame)
    if not path.is_dir():
        return Annotation.from_long_table(pd.DataFrame())
    if dir_looks_like_annotation(path):
        return Annotation.from_annotation_dir(path)
    tables: Dict[str, pd.DataFrame] = {}
    for csv in sorted(path.glob("*.csv")):
        if csv.name.startswith("."):
            continue
        try:
            frame = pd.read_csv(csv)
        except Exception:
            continue
        if is_abundance_table(frame):
            tables[csv.stem] = frame
    if tables:
        return Annotation.from_abundance_tables(tables)
    return Annotation.from_annotation_dir(path)


def observed_n_samples(data: Any) -> int:
    tables = input_to_abundance_tables(data)
    for table in tables.values():
        n = n_sample_columns(table)
        if n:
            return max(len(n), 1)
    frame = getattr(data, "DataFrame", data)
    if isinstance(frame, pd.DataFrame) and "sample" in frame.columns and len(frame.index):
        return max(int(pd.Series(frame["sample"].astype(str)).nunique()), 1)
    return 1


def skip_abundance_filename(name: str) -> bool:
    low = str(name).lower()
    if low.startswith("."):
        return True
    if low == "table_selection.json":
        return True
    if low.startswith("combined_annotation_table"):
        return True
    return False


def abundance_csv_paths(directory: PathLike) -> List[Path]:
    folder = Path(directory)
    if not folder.is_dir():
        return []
    return [
        path
        for path in sorted(folder.glob("*.csv"))
        if path.is_file() and not skip_abundance_filename(path.name)
    ]


def load_abundance_dir(directory: PathLike) -> Dict[str, pd.DataFrame]:
    tables: Dict[str, pd.DataFrame] = {}
    for path in abundance_csv_paths(directory):
        try:
            frame = pd.read_csv(path)
        except Exception:
            continue
        if is_abundance_table(frame):
            tables[path.stem] = normalize_abundance_table(frame)
    return tables


def has_abundance_tables(directory: PathLike) -> bool:
    return bool(load_abundance_dir(directory))


def observed_abundance_dir(output_dir: PathLike) -> Path:
    return Path(output_dir) / INITIAL_ABUNDANCE_DIR


def regenerated_abundance_dir(output_dir: PathLike) -> Path:
    return Path(output_dir) / REGENERATED_ABUNDANCE_DIR


def write_abundance_dir(dest: PathLike, tables: Dict[str, pd.DataFrame]) -> Path:
    from samovar.regenerate import _write_abundance_tables

    out = Path(dest)
    _write_abundance_tables(out, tables)
    return out


def collect_observed_abundance(output_dir: PathLike) -> Dict[str, pd.DataFrame]:
    """Find observed abundance in the run dir (already staged, annotations, or loose CSVs)."""
    root = Path(output_dir)
    staged = load_abundance_dir(observed_abundance_dir(root))
    if staged:
        return staged
    for folder in (root / "abundance", root / "initial_annotations", root):
        tables = load_abundance_dir(folder)
        if tables:
            return tables
    ann = root / "initial_annotations"
    if ann.is_dir():
        loaded = load_table_input(ann)
        tables = input_to_abundance_tables(loaded)
        if tables:
            return tables
    return {}


def materialize_observed_abundance(output_dir: PathLike) -> Dict[str, pd.DataFrame]:
    """Write canonical ``initial_abundance/*.csv`` from annotations or OTU tables."""
    tables = collect_observed_abundance(output_dir)
    if not tables:
        raise FileNotFoundError(
            f"No abundance or annotation tables under {output_dir}. "
            "Put OTU/abundance CSVs in initial_abundance/ or annotation CSVs in initial_annotations/."
        )
    write_abundance_dir(observed_abundance_dir(output_dir), tables)
    return tables


def convert_to_abundance_dir(
    source: PathLike,
    dest: PathLike,
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, pd.DataFrame]:
    """Write annotator CSVs (``taxid`` + ``N_<sample>``) from annotations or OTU tables.

    If ``dest`` already has abundance CSVs, they are reused. Generative table
    regenerators run only when ``config`` asks for a mode other than ``direct``.
    """
    dest_p = Path(dest)
    if has_abundance_tables(dest_p):
        return load_abundance_dir(dest_p)
    src = Path(source)
    if not src.exists():
        raise FileNotFoundError(f"No annotation or abundance source at {source}")
    cfg = dict(config or {})
    from samovar.regenerate import regenerate_annotation_tables, resolve_regeneration_mode

    raw_mode = (
        cfg.get("table_reads_generators")
        or cfg.get("table_reads_generator")
        or cfg.get("regeneration_mode")
        or "direct"
    )
    if isinstance(raw_mode, (list, tuple)):
        first = raw_mode[0] if raw_mode else "direct"
    else:
        first = raw_mode
    kind, canon = resolve_regeneration_mode(first)
    generative = not (kind == "builtin" and canon == "direct") or (
        isinstance(raw_mode, (list, tuple)) and len(raw_mode) > 1
    )
    if generative:
        return regenerate_annotation_tables(src, dest_p, cfg, select_best=False)
    tables = load_abundance_dir(src) if src.is_dir() else {}
    if not tables:
        loaded = load_table_input(src)
        tables = input_to_abundance_tables(loaded)
    if not tables:
        raise FileNotFoundError(
            f"No abundance or annotation tables under {source}. "
            "Need OTU CSVs (taxid + N_*) or long *.annotation.csv with taxID_*."
        )
    write_abundance_dir(dest_p, tables)
    return tables


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(prog="python -m samovar.abundance")
    sub = parser.add_subparsers(dest="command", required=True)
    mat = sub.add_parser("materialize", help="Write outdir/initial_abundance from annotations or OTU CSVs")
    mat.add_argument("--output_dir", "--outdir", dest="output_dir", required=True)
    conv = sub.add_parser(
        "convert",
        help="annotation/OTU directory → abundance CSVs (annotation2abundance)",
    )
    conv.add_argument("--source", "--annotation_dir", dest="source", required=True)
    conv.add_argument("--dest", "--abundance_dir", dest="dest", required=True)
    conv.add_argument("--config", default="", help="Optional YAML (table regenerator settings)")
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.command == "materialize":
        tables = materialize_observed_abundance(args.output_dir)
        print(f"Wrote {len(tables)} abundance table(s) under {observed_abundance_dir(args.output_dir)}")
        return 0
    if args.command == "convert":
        cfg: Dict[str, Any] = {}
        if args.config:
            import yaml

            cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
        tables = convert_to_abundance_dir(args.source, args.dest, cfg)
        print(f"Wrote {len(tables)} abundance table(s) under {args.dest}")
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
