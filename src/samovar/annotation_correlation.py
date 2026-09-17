"""Spearman correlation of tax / feat / true columns on an Annotation table.

This is a pipeline helper (not a ``samovar tools import`` contract). It replaces
the old pairwise annotator cross-tab heatmaps with one generalised matrix.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Union

import numpy as np
import pandas as pd

from samovar.annotation_columns import feat_annotator_columns, is_feat_column, is_tax_column

PathLike = Union[str, Path]
MAX_FEAT_COLUMNS = 20
CORRELATION_TYPES = frozenset(
    {"cv", "cross-validation", "correlation", "spearman", "corr"}
)
DIVERGING_COLSTOPS = [
    [0.0, "#2166AC"],
    [0.25, "#92C5DE"],
    [0.5, "#F7F7F7"],
    [0.75, "#F4A582"],
    [1.0, "#B2182B"],
]


def annotation_stage_allows_correlation(
    annotation_dir: Optional[PathLike] = None,
    output_dir: Optional[PathLike] = None,
) -> bool:
    """True for regenerated / reprofiled viz folders (not Raw / initial)."""
    blob = f"{annotation_dir or ''} {output_dir or ''}".lower().replace("-", "_")
    return "regenerated" in blob or "reprofiled" in blob


def _as_frame(annotation: Any) -> pd.DataFrame:
    if annotation is None:
        return pd.DataFrame()
    if isinstance(annotation, pd.DataFrame):
        return annotation
    to_long = getattr(annotation, "to_long_table", None) or getattr(
        annotation, "to_dataframe", None
    )
    if callable(to_long):
        frame = to_long()
        if isinstance(frame, pd.DataFrame):
            return frame
    table = getattr(annotation, "table", None)
    if isinstance(table, pd.DataFrame):
        return table
    if isinstance(annotation, (str, Path)) and Path(annotation).is_file():
        path = Path(annotation)
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        return pd.read_csv(path, sep=sep)
    return pd.DataFrame()


def _spearman_pair(left: pd.Series, right: pd.Series) -> float:
    """Spearman ρ; 0 when ranks cannot be compared (constant or too few pairs)."""
    both = pd.concat([left, right], axis=1, keys=["a", "b"]).dropna()
    if len(both) < 2:
        return 0.0
    if int(both["a"].nunique()) < 2 or int(both["b"].nunique()) < 2:
        return 0.0
    rho = both["a"].corr(both["b"], method="spearman")
    if rho is None or pd.isna(rho):
        return 0.0
    return float(rho)


def _rankable_series(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    if float(numeric.notna().mean()) >= 0.5:
        return numeric
    extracted = series.astype(str).str.extract(r"(-?\d+)", expand=False)
    from_digits = pd.to_numeric(extracted, errors="coerce")
    if from_digits.notna().any():
        return from_digits
    text = series.astype(str).str.strip()
    blank = text.str.lower().isin({"", "nan", "none", "na", "<na>"})
    codes, _uniques = pd.factorize(text.where(~blank), sort=True)
    out = pd.Series(codes, index=series.index, dtype=float)
    valid = pd.Series(np.asarray(codes) >= 0, index=series.index)
    return out.where(valid)


META_SKIP = frozenset(
    {"seq", "sample", "true", "read_type", "length", "classified", "taxa"}
)


def correlation_columns(annotation: Any) -> List[str]:
    frame = _as_frame(annotation)
    tax: List[str] = []
    for col in frame.columns:
        if is_tax_column(str(col), include_true=False):
            tax.append(str(col))
    if not tax:
        for col in frame.columns:
            name = str(col)
            low = name.lower()
            if low in META_SKIP or "conf" in low:
                continue
            if is_feat_column(name):
                continue
            tax.append(name)
    feat = [str(c) for c in feat_annotator_columns(frame.columns)]
    cols = tax + feat
    if "true" in frame.columns and "true" not in cols:
        cols.append("true")
    seen = set()
    ordered = []
    for col in cols:
        if col not in seen and col in {str(c) for c in frame.columns}:
            seen.add(col)
            ordered.append(col)
    return ordered


def importance_key(name: str) -> str:
    """Match FI names to heatmap columns (taxid_/feat_ prefixes and run indices)."""
    text = str(name).strip().lower().replace("-", "_")
    text = re.sub(r"^(taxid|tax|feat)_+", "", text)
    text = re.sub(r"^([a-z][a-z0-9]*)_\d+_", r"\1_", text)
    text = re.sub(r"^([a-z][a-z0-9]*)_\d+$", r"\1", text)
    return text


def load_feature_importance_scores(root: Optional[PathLike]) -> Optional[pd.Series]:
    """Load the default model's feature-importance vector from a run directory."""
    if not root:
        return None
    base = Path(root)
    folders: List[Path] = []
    low = base.name.lower()
    if low.endswith("_plots") or low.endswith("_annotations"):
        folders.append(base.parent / "feature_importance_plots")
        folders.append(base.parent / "reprofiled_annotations" / "feature_importance")
        folders.append(base.parent)
    folders.extend(
        [
            base / "feature_importance_plots",
            base / "reprofiled_annotations" / "feature_importance",
            base,
        ]
    )
    files: List[Path] = []
    seen = set()
    for folder in folders:
        key = str(folder)
        if key in seen:
            continue
        seen.add(key)
        if folder.is_dir():
            files.extend(sorted(p for p in folder.glob("feature_importance_*.tsv") if p.is_file()))
    preferred = [
        path
        for path in files
        if any(tok in path.stem.lower() for tok in ("best", "randomforest", "rf"))
    ]
    for path in preferred + [p for p in files if p not in preferred]:
        try:
            table = pd.read_csv(path, sep="\t")
        except Exception:
            continue
        if "feature" not in table.columns or "importance" not in table.columns:
            continue
        series = pd.Series(
            pd.to_numeric(table["importance"], errors="coerce").fillna(0.0).to_numpy(dtype=float),
            index=table["feature"].astype(str),
        )
        if not series.empty:
            return series
    return None


def _importance_lookup(scores: Optional[pd.Series]) -> Dict[str, float]:
    lookup: Dict[str, float] = {}
    if scores is None:
        return lookup
    for name, value in scores.items():
        key = importance_key(str(name))
        if not key:
            continue
        num = 0.0 if pd.isna(value) else float(value)
        lookup[key] = max(lookup.get(key, 0.0), num)
    return lookup


def _column_rank_score(
    name: str,
    lookup: Mapping[str, float],
    rho_abs: Mapping[str, float],
) -> float:
    if lookup:
        return float(lookup.get(importance_key(name), 0.0))
    return float(rho_abs.get(name, 0.0))


def _select_and_order(
    columns: Sequence[str],
    numeric: pd.DataFrame,
    max_feat: int,
    importance: Optional[pd.Series] = None,
) -> List[str]:
    feat_set = {
        c
        for c in columns
        if c in feat_annotator_columns(columns) or str(c).lower().startswith("feat_")
    }
    true_cols = [c for c in columns if str(c).lower() == "true"]
    core = [c for c in columns if c not in feat_set and str(c).lower() != "true"]
    feat = [c for c in columns if c in feat_set]
    rho_abs: Dict[str, float] = {}
    if "true" in numeric.columns:
        target = numeric["true"]
        for col in list(core) + feat:
            rho_abs[col] = abs(_spearman_pair(numeric[col], target))
    else:
        for col in feat:
            rho_abs[col] = float(numeric[col].var(skipna=True) or 0.0)
    lookup = _importance_lookup(importance)
    feat_ranked = sorted(
        feat,
        key=lambda col: (_column_rank_score(col, lookup, rho_abs), str(col)),
        reverse=True,
    )
    keep_feat = feat_ranked[: max(0, int(max_feat))]
    selected = core + keep_feat
    selected = sorted(
        selected,
        key=lambda col: (_column_rank_score(col, lookup, rho_abs), str(col)),
        reverse=True,
    )
    return selected + true_cols


def spearman_annotation_correlation(
    annotation,
    *,
    max_feat: int = MAX_FEAT_COLUMNS,
    importance: Optional[pd.Series] = None,
    importance_dir: Optional[PathLike] = None,
) -> pd.DataFrame:
    """Pairwise Spearman ρ among tax annotators, top feat_ columns, and true taxid."""
    frame = _as_frame(annotation)
    names = correlation_columns(frame)
    if len(names) < 2:
        return pd.DataFrame()
    numeric = pd.DataFrame({name: _rankable_series(frame[name]) for name in names})
    scores = importance if importance is not None else load_feature_importance_scores(importance_dir)
    names = _select_and_order(names, numeric, max_feat, scores)
    if len(names) < 2:
        return pd.DataFrame()
    numeric = numeric.loc[:, names]
    n = len(names)
    arr = np.eye(n, dtype=float)
    for i, left in enumerate(names):
        for j in range(i + 1, n):
            rho = _spearman_pair(numeric[left], numeric[names[j]])
            arr[i, j] = arr[j, i] = rho
    return pd.DataFrame(arr, index=names, columns=names)


def _save_cns_heatmap(matrix: pd.DataFrame, path: Path, title: str) -> None:
    try:
        import cnsplots as cns
        import matplotlib.pyplot as plt
        from matplotlib.colors import TwoSlopeNorm

        fig, ax = plt.subplots(
            figsize=(max(4.0, 0.38 * matrix.shape[1] + 2.2), max(4.0, 0.38 * matrix.shape[0] + 2.2))
        )
        values = matrix.to_numpy(dtype=float)
        norm = TwoSlopeNorm(vmin=-1.0, vcenter=0.0, vmax=1.0)
        im = ax.imshow(values, cmap="RdBu_r", norm=norm, origin="upper", aspect="auto")
        ax.set_xticks(range(len(matrix.columns)))
        ax.set_yticks(range(len(matrix.index)))
        ax.set_xticklabels([str(c) for c in matrix.columns], rotation=90, fontsize=7)
        ax.set_yticklabels([str(c) for c in matrix.index], fontsize=7)
        ax.set_title(title)
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="Spearman ρ")
        n = min(matrix.shape[0], 40)
        if matrix.size <= 1600:
            for i in range(n):
                for j in range(min(matrix.shape[1], 40)):
                    val = values[i, j]
                    if val == val:
                        ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=5)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        save = getattr(cns, "savefig", None)
        if callable(save):
            save(str(path))
            try:
                save(str(path.with_suffix(".svg")))
            except Exception:
                fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
        else:
            fig.savefig(path, dpi=150, bbox_inches="tight")
            fig.savefig(path.with_suffix(".svg"), bbox_inches="tight")
        plt.close(fig)
    except Exception:
        pass


def _save_altair_heatmap(matrix: pd.DataFrame, path: Path, title: str) -> None:
    try:
        import altair as alt

        long = (
            matrix.rename_axis("y")
            .reset_index()
            .melt(id_vars="y", var_name="x", value_name="rho")
        )
        chart = (
            alt.Chart(long)
            .mark_rect()
            .encode(
                x=alt.X("x:N", title="Column", sort=list(matrix.columns.astype(str))),
                y=alt.Y("y:N", title="Column", sort=list(matrix.index.astype(str))),
                color=alt.Color(
                    "rho:Q",
                    scale=alt.Scale(domain=[-1, 0, 1], scheme="redblue", reverse=True),
                    title="Spearman ρ",
                ),
                tooltip=["x", "y", "rho"],
            )
            .properties(title=title, width=max(180, 18 * matrix.shape[1]), height=max(180, 18 * matrix.shape[0]))
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        chart.save(str(path))
    except Exception:
        pass


def write_spearman_correlation(
    matrix: pd.DataFrame,
    dest: PathLike,
    *,
    stem: str = "spearman_correlation",
    section_name: str = "Spearman correlation",
    description: str = "",
) -> Dict[str, Any]:
    """Write MultiQC heatmap JSON plus optional altair / cnsplots figures."""
    from samovar.stage_report import write_heatmap_mqc

    folder = Path(dest)
    folder.mkdir(parents=True, exist_ok=True)
    if matrix is None or matrix.empty:
        return {}
    tsv = folder / f"{stem}.tsv"
    matrix.to_csv(tsv, sep="\t")
    desc = description or (
        "Spearman rank correlation of tax annotator columns, the top 20 feat_ "
        "columns by feature importance (default), and true taxid (when present)."
    )
    write_heatmap_mqc(
        matrix,
        folder / f"{stem}_mqc.json",
        section_name=section_name,
        description=desc,
        xlab="Column",
        ylab="Column",
        min_value=-1,
        max_value=1,
        colstops=DIVERGING_COLSTOPS,
    )
    _save_altair_heatmap(matrix, folder / f"{stem}.html", section_name)
    _save_cns_heatmap(matrix, folder / f"{stem}.png", section_name)
    return {"matrix": matrix, "path": str(folder / f"{stem}_mqc.json")}


def plot_spearman_annotation_correlation(
    annotation,
    output_dir: PathLike,
    *,
    read_type: Optional[str] = None,
    max_feat: int = MAX_FEAT_COLUMNS,
    importance: Optional[pd.Series] = None,
    importance_dir: Optional[PathLike] = None,
) -> Optional[pd.DataFrame]:
    matrix = spearman_annotation_correlation(
        annotation,
        max_feat=max_feat,
        importance=importance,
        importance_dir=importance_dir if importance_dir is not None else output_dir,
    )
    if matrix.empty:
        return None
    stem = "spearman_correlation"
    title = "Spearman correlation (tax / top-20 feat by importance / true)"
    if read_type:
        stem = f"{stem}.{read_type}"
        title = f"{title} — {read_type}"
    write_spearman_correlation(matrix, output_dir, stem=stem, section_name=title)
    return matrix


def refresh_spearman_heatmaps(run_dir: PathLike) -> None:
    """Rewrite Regenerated/Reprofiled Spearman heatmaps using saved feature importance."""
    from samovar.annotation_io import read_annotation_dir
    from samovar.viz_annotation import _strip_annotator_prefixes

    root = Path(run_dir)
    scores = load_feature_importance_scores(root)
    for stage in ("regenerated", "reprofiled"):
        annotation_dir = root / f"{stage}_annotations"
        dest = root / f"{stage}_annotations_plots"
        if not annotation_dir.is_dir():
            continue
        try:
            data = read_annotation_dir(str(annotation_dir))
        except Exception:
            data = pd.DataFrame()
        if data.empty:
            combined = annotation_dir / "combined_annotation_table.csv"
            if combined.is_file():
                data = pd.read_csv(combined)
        if data.empty:
            continue
        work = _strip_annotator_prefixes(data)
        plot_spearman_annotation_correlation(work, dest, importance=scores)
