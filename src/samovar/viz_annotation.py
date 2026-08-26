"""Python analog of ``R/viz_annotation.R``.

Publication PNGs use cnsplots (+ matplotlib). Interactive HTML uses Altair.
Rank remapping is in-memory only; combined CSVs on disk stay at exact taxIDs.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from samovar.annotation_io import annotator_columns, read_annotation_dir
from samovar.parse_annotators import (
    ensure_taxid_name_map,
    remap_taxid_dataframe,
)

PALETTE_F1 = [
    "#FFFFFF",
    "#FFFFE5",
    "#F7FCB9",
    "#E6E487",
    "#D9F0A3",
    "#ADDD8E",
    "#78C679",
    "#41AB5D",
    "#238443",
    "#006837",
    "#004529",
]


def _use_agg_backend() -> None:
    import matplotlib

    matplotlib.use("Agg")


def require_viz_backend() -> str:
    """Return ``cnsplots`` (preferred) or ``altair``.

    Pipeline plots need at least one of these. ``./install.sh`` installs both
    via pyproject; a compute node without either must fail rather than skip.
    """
    try:
        import cnsplots  # noqa: F401

        return "cnsplots"
    except ImportError:
        pass
    try:
        import altair  # noqa: F401

        return "altair"
    except ImportError:
        pass
    raise RuntimeError(
        "Visualization requires cnsplots (preferred) or altair. "
        "Reinstall with ./install.sh or: pip install 'cnsplots>=0.6.0' altair"
    )


def _setup_cns():
    _use_agg_backend()
    try:
        import cnsplots as cns

        if hasattr(cns, "setup_matplotlib"):
            cns.setup_matplotlib()
        return cns
    except Exception:
        return None


def is_special_taxon(x) -> bool:
    token = str(x).strip().lower()
    return token in {"0", "other", "unclassified", "none", "na", "nan", ""}


def normalize_taxon_token(x) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "0"
    text = str(x).strip()
    if text.lower() in {"", "nan", "none", "na", "other"}:
        return "other" if text.lower() == "other" else "0"
    if text.endswith(".0") and text[:-2].lstrip("-").isdigit():
        text = text[:-2]
    return text


def fpc_taxon_order(freq_df: pd.DataFrame, col_x: str, col_y: str, freq_col: str = "Freq") -> List[str]:
    x = freq_df[col_x].map(normalize_taxon_token)
    y = freq_df[col_y].map(normalize_taxon_token)
    f = pd.to_numeric(freq_df[freq_col], errors="coerce").fillna(0)
    core = [t for t in pd.unique(pd.concat([x, y], ignore_index=True)) if not is_special_taxon(t)]
    if len(core) <= 1:
        return core
    idx = {t: i for i, t in enumerate(core)}
    mat = np.zeros((len(core), len(core)), dtype=float)
    for a, b, n in zip(y, x, f):
        if a in idx and b in idx:
            mat[idx[a], idx[b]] += float(n)
    mat_s = mat + mat.T
    np.fill_diagonal(mat_s, np.diag(mat))
    try:
        if mat_s.shape[0] < 2 or np.std(mat_s) == 0:
            order = np.argsort(-mat_s.sum(axis=1))
        else:
            centered = mat_s - mat_s.mean(axis=0)
            _, _, vt = np.linalg.svd(centered, full_matrices=False)
            order = np.argsort(centered @ vt[0])
    except Exception:
        order = np.argsort(-mat_s.sum(axis=1))
    return [core[i] for i in order]


def axis_levels_from_fpc(present: Iterable, fpc_core: Sequence[str]) -> List[str]:
    present = [normalize_taxon_token(v) for v in pd.Series(list(present), dtype=str).unique()]
    core = [t for t in fpc_core if t in present]
    extra = [t for t in present if not is_special_taxon(t) and t not in core]
    spec = [t for t in ("0", "other") if t in present]
    return core + extra + spec


def _strip_annotator_prefixes(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    rename = {}
    for col in out.columns:
        name = re.sub(r"^taxid_", "", str(col), flags=re.I)
        name = re.sub(r"^N_", "", name, flags=re.I)
        rename[col] = name
    out = out.rename(columns=rename)
    if "length" not in out.columns:
        out["length"] = np.nan
    if "sample" not in out.columns:
        out["sample"] = np.nan
    return out


def _annotator_names(df: pd.DataFrame) -> List[str]:
    skip = {"length", "sample", "true", "seq", "read_type"}
    names = []
    for col in df.columns:
        if col in skip or "confidence" in str(col).lower() or str(col).endswith("_conf"):
            continue
        if str(col).startswith("taxID") or str(col) in skip:
            continue
        names.append(col)
    # After prefix strip, annotators are remaining non-meta columns that were taxID_/N_
    if names:
        return names
    return [c for c in df.columns if c not in skip and "conf" not in str(c).lower()]


def _collapse_other(pred: pd.Series, true: pd.Series) -> pd.Series:
    allowed = set(true.map(normalize_taxon_token)) | {"0"}
    out = pred.map(normalize_taxon_token)
    return out.where(out.isin(allowed), "other")


def _confusion_counts(true: pd.Series, pred: pd.Series) -> pd.DataFrame:
    tmp = pd.DataFrame(
        {
            "true": true.map(normalize_taxon_token),
            "pred": pred.map(normalize_taxon_token),
        }
    )
    return tmp.value_counts().rename("Freq").reset_index()


def _f1_from_counts(counts: pd.DataFrame) -> float:
    rect = counts["true"].astype(str) == counts["pred"].astype(str)
    tp = counts.loc[rect, "Freq"].sum()
    off = counts.loc[~rect, "Freq"].sum()
    precision = tp / (tp + off) if (tp + off) else 0.0
    recall = precision  # same off-diagonal pool as the R implementation
    if precision + recall == 0:
        return 0.0
    return float(2 * precision * recall / (precision + recall))


def _r2_table(true: pd.Series, pred: pd.Series) -> Tuple[pd.DataFrame, float]:
    counts = _confusion_counts(true, pred)
    true_tot = counts.groupby("true", as_index=False)["Freq"].sum().rename(columns={"Freq": "true_n"})
    pred_tot = counts.groupby("pred", as_index=False)["Freq"].sum().rename(columns={"Freq": "pred_n"})
    pred_tot = pred_tot.rename(columns={"pred": "true"})
    table = true_tot.merge(pred_tot, on="true", how="left").fillna(0)
    if table["pred_n"].sum() <= 0:
        return table, float("nan")
    ss_tot = ((table["true_n"] - table["true_n"].mean()) ** 2).sum()
    ss_res = ((table["pred_n"] - table["true_n"]) ** 2).sum()
    r2 = 0.0 if ss_tot == 0 else float(1 - ss_res / ss_tot)
    return table, r2


def _clip_taxon_label(text: str, max_len: int = 25) -> str:
    s = str(text)
    if len(s) <= max_len:
        return s
    return s[: max_len - 1] + "…"


def _taxon_label(token: str, name_map: Dict[str, str]) -> str:
    if is_special_taxon(token):
        raw = "other" if str(token).lower() == "other" else "0"
        return _clip_taxon_label(raw)
    key = str(token)
    name = name_map.get(key, name_map.get(key.rstrip(".0"), key))
    return _clip_taxon_label(name)


def _trim_matrix(matrix: pd.DataFrame, max_axis: int = 40) -> pd.DataFrame:
    """Keep abundant taxa so heatmaps stay publication-sized (avoids OOM)."""
    if matrix.shape[0] <= max_axis and matrix.shape[1] <= max_axis:
        return matrix
    row_keep = list(matrix.sum(axis=1).sort_values(ascending=False).head(max_axis).index)
    col_keep = list(matrix.sum(axis=0).sort_values(ascending=False).head(max_axis).index)
    special = [t for t in ("0", "other") if t in matrix.index or t in matrix.columns]
    rows = []
    for t in row_keep + special:
        if t in matrix.index and t not in rows:
            rows.append(t)
    cols = []
    for t in col_keep + special:
        if t in matrix.columns and t not in cols:
            cols.append(t)
    return matrix.reindex(index=rows[:max_axis], columns=cols[:max_axis], fill_value=0)


def _align_xticklabels(ax, labels: Sequence[str], italic: bool) -> None:
    """Draw 90° x labels flush with the spine (ggplot hjust=1, vjust=0.5).

    Native tick labels + cnsplots ``xtick.alignment`` keep the string centered on
    the axis, so half of each name crosses into the heatmap. Explicit ``ax.text``
    at y=0 in axes coordinates puts the top of the text on the spine.
    """
    ax.set_xticklabels([])
    ax.tick_params(axis="x", pad=2, length=3, width=0.6, labelbottom=False)
    for artist in list(ax.texts):
        if getattr(artist, "_samovar_xtick", False):
            artist.remove()
    from matplotlib.transforms import ScaledTranslation

    # A few points below the spine so glyphs do not cross into the heatmap.
    offset = ScaledTranslation(0, -8 / 72.0, ax.figure.dpi_scale_trans)
    trans = ax.get_xaxis_transform() + offset
    for j, text in enumerate(labels):
        t = ax.text(
            j,
            0.0,
            text,
            transform=trans,
            rotation=90,
            ha="center",
            va="top",
            fontsize=7,
            clip_on=False,
            zorder=6,
        )
        t._samovar_xtick = True
        if italic and text not in {"0", "other"}:
            t.set_style("italic")


def _xtick_text_artists(ax):
    return [t for t in ax.texts if getattr(t, "_samovar_xtick", False)]


def _outline_diagonal_cells(ax, matrix: pd.DataFrame) -> None:
    from matplotlib.patches import Rectangle

    for i, y_tok in enumerate(matrix.index):
        for j, x_tok in enumerate(matrix.columns):
            if str(y_tok) != str(x_tok):
                continue
            ax.add_patch(
                Rectangle(
                    (j - 0.5, i - 0.5),
                    1,
                    1,
                    fill=False,
                    edgecolor="0.45",
                    linewidth=2.0,
                    zorder=4,
                    clip_on=False,
                )
            )


def _annotate_heatmap_counts(ax, matrix: pd.DataFrame, log_values: np.ndarray) -> None:
    if max(matrix.shape) > 16:
        return
    vmin = float(np.nanmin(log_values)) if log_values.size else 0.0
    vmax = float(np.nanmax(log_values)) if log_values.size else 1.0
    span = vmax - vmin
    n_pal = max(2, len(PALETTE_F1))
    dark_cut = (n_pal - 3) / (n_pal - 1)
    for i, y_tok in enumerate(matrix.index):
        for j, x_tok in enumerate(matrix.columns):
            val = matrix.iat[i, j]
            if val <= 0:
                continue
            frac = 0.0 if span <= 0 else (log_values[i, j] - vmin) / span
            unclassified = is_special_taxon(x_tok) or is_special_taxon(y_tok)
            if frac >= dark_cut:
                color = "white"
            elif unclassified:
                color = "saddlebrown"
            else:
                color = "black"
            ax.text(j, i, str(int(val)), ha="center", va="center", fontsize=6, color=color, zorder=5)


def _save_heatmap_png(
    matrix: pd.DataFrame,
    path: Path,
    title: str,
    xlabel: str,
    ylabel: str,
    caption: str,
    name_map: Dict[str, str],
    italic: bool,
) -> None:
    cns = _setup_cns()
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    x_labels = [_taxon_label(c, name_map) for c in matrix.columns]
    y_labels = [_taxon_label(c, name_map) for c in matrix.index]
    values = np.log10(matrix.to_numpy(dtype=float) + 1)
    n = max(matrix.shape)
    longest_x = max((len(s) for s in x_labels), default=1)
    longest_y = max((len(s) for s in y_labels), default=1)
    # Extra height/width so rotated x ticks and y names are not clipped.
    width = max(5.0, 0.42 * n + 1.4 + 0.07 * min(longest_y, 24))
    height = max(5.0, 0.42 * n + 1.8 + 0.08 * min(longest_x, 25))
    fig, ax = plt.subplots(figsize=(width, height), facecolor="white")
    ax.set_facecolor("white")
    ax.margins(0)
    cmap = LinearSegmentedColormap.from_list("samovar_f1", PALETTE_F1)
    im = ax.imshow(values, cmap=cmap, aspect="equal", origin="upper")
    ax.set_xlim(-0.5, max(len(x_labels) - 0.5, 0.5))
    ax.set_ylim(max(len(y_labels) - 0.5, 0.5), -0.5)
    ax.set_xticks(range(len(x_labels)))
    ax.set_yticks(range(len(y_labels)))
    ax.set_xticklabels([])
    ax.set_yticklabels(y_labels, fontsize=7)
    if italic:
        for lab in ax.get_yticklabels():
            if lab.get_text() not in {"0", "other"}:
                lab.set_style("italic")
    _align_xticklabels(ax, x_labels, italic)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel, labelpad=8)
    ax.set_title(title)
    ax.set_xticks(np.arange(-0.5, len(x_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(y_labels), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.3)
    ax.tick_params(which="minor", bottom=False, left=False)
    _outline_diagonal_cells(ax, matrix)
    _annotate_heatmap_counts(ax, matrix, values)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    dpi = float(fig.dpi)
    tick_artists = _xtick_text_artists(ax)
    tick_h = max(
        (t.get_window_extent(renderer=renderer).height for t in tick_artists),
        default=0.0,
    )
    ax.xaxis.labelpad = tick_h * 72.0 / dpi + 10
    fig.tight_layout()
    _align_xticklabels(ax, x_labels, italic)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    if caption:
        boxes = [ax.xaxis.label.get_window_extent(renderer=renderer)]
        boxes.extend(t.get_window_extent(renderer=renderer) for t in _xtick_text_artists(ax))
        lowest = min(b.y0 for b in boxes)
        _, y_fig = fig.transFigure.inverted().transform((0, lowest))
        fig.text(
            0.08,
            y_fig - 0.08,
            caption,
            ha="left",
            va="top",
            fontsize=8,
            color="black",
            transform=fig.transFigure,
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path,
        dpi=150,
        bbox_inches="tight",
        pad_inches=0.5,
        facecolor="white",
        transparent=False,
        edgecolor="none",
    )
    plt.close("all")


def _save_scatter_png(table: pd.DataFrame, path: Path, title: str, caption: str) -> None:
    cns = _setup_cns()
    import matplotlib.pyplot as plt

    if cns is not None and hasattr(cns, "figure"):
        cns.figure(220, 220)
        ax = plt.gca()
    else:
        _, ax = plt.subplots(figsize=(8, 8))
    x = table["pred_n"].to_numpy(dtype=float)
    y = table["true_n"].to_numpy(dtype=float)
    ratio = np.log10((y + 1e-9) / (x + 1e-9))
    ax.scatter(x, y, c=ratio, cmap="RdYlGn_r", edgecolors="none")
    if len(x) >= 2:
        coef = np.polyfit(x, y, 1)
        xs = np.linspace(min(x.min(), y.min()), max(x.max(), y.max()), 50)
        ax.plot(xs, np.polyval(coef, xs), color="0.4", alpha=0.5)
    lims = [0, max(x.max(), y.max(), 1) * 1.05]
    ax.plot(lims, lims, color="red", linestyle="--")
    for _, row in table.iterrows():
        ax.text(row["pred_n"], row["true_n"] + max(y.max() / 20, 1), str(row["true"]), fontsize=7)
    ax.set_xlabel("Predicted taxon", labelpad=8)
    ax.set_ylabel("True taxon")
    ax.set_title(title)
    ax.set_aspect("equal", adjustable="box")
    fig = ax.figure
    fig.tight_layout()
    if caption:
        ax.annotate(
            caption,
            xy=(0.0, 0.0),
            xycoords="axes fraction",
            xytext=(0, -56),
            textcoords="offset points",
            ha="left",
            va="top",
            fontsize=8,
            annotation_clip=False,
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        path,
        dpi=150,
        bbox_inches="tight",
        pad_inches=0.35,
        facecolor="white",
        transparent=False,
        edgecolor="none",
    )
    plt.close("all")


def _save_confidence_png(df: pd.DataFrame, path: Path, title: str) -> None:
    cns = _setup_cns()
    import matplotlib.pyplot as plt

    if cns is not None and hasattr(cns, "figure"):
        cns.figure(240, 180)
        ax = plt.gca()
    else:
        _, ax = plt.subplots(figsize=(8, 6))
    taxa = list(df["y"].astype(str).unique())
    data = [df.loc[df["y"].astype(str) == t, "x"].dropna().to_numpy() for t in taxa]
    try:
        ax.boxplot(data, vert=False, tick_labels=taxa, showfliers=False)
    except TypeError:
        ax.boxplot(data, vert=False, labels=taxa, showfliers=False)
    rng = np.random.default_rng(0)
    for i, t in enumerate(taxa, start=1):
        xs = df.loc[df["y"].astype(str) == t, "x"].dropna().to_numpy()
        ax.scatter(xs, np.full_like(xs, i, dtype=float) + rng.uniform(-0.08, 0.08, size=len(xs)), s=4, alpha=0.15)
    ax.set_xlabel("confidence")
    ax.set_ylabel("taxID")
    ax.set_title(title)
    path.parent.mkdir(parents=True, exist_ok=True)
    if cns is not None and hasattr(cns, "savefig"):
        cns.savefig(str(path))
    else:
        plt.tight_layout()
        plt.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")


def _altair_heatmap(matrix: pd.DataFrame, title: str, xlabel: str, ylabel: str):
    import altair as alt

    long = matrix.reset_index().melt(id_vars=matrix.index.name or "index", var_name="x", value_name="Freq")
    ycol = long.columns[0]
    long = long.rename(columns={ycol: "y"})
    return (
        alt.Chart(long)
        .mark_rect()
        .encode(
            x=alt.X("x:N", title=xlabel),
            y=alt.Y("y:N", title=ylabel, sort="descending"),
            color=alt.Color("Freq:Q", scale=alt.Scale(scheme="yellowgreen")),
            tooltip=["x", "y", "Freq"],
        )
        .properties(title=title)
    )


def _altair_scatter(table: pd.DataFrame, title: str):
    import altair as alt

    return (
        alt.Chart(table)
        .mark_point()
        .encode(
            x=alt.X("pred_n:Q", title="Predicted taxon"),
            y=alt.Y("true_n:Q", title="True taxon"),
            tooltip=["true", "true_n", "pred_n"],
        )
        .properties(title=title)
    )


def _save_altair(chart, path: Path, n_cells: int = 0) -> None:
    if n_cells > 2500:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        chart.save(str(path))
    except Exception:
        try:
            path.with_suffix(".json").write_text(chart.to_json())
        except Exception:
            pass


def _pivot_matrix(counts: pd.DataFrame, x_col: str, y_col: str, x_levels: List[str], y_levels: List[str]) -> pd.DataFrame:
    mat = (
        counts.pivot_table(index=y_col, columns=x_col, values="Freq", aggfunc="sum", fill_value=0)
        .reindex(index=y_levels, columns=x_levels, fill_value=0)
    )
    mat.index.name = y_col
    return mat


def viz_annotation(
    data: pd.DataFrame,
    type: Sequence[str] = ("f1", "R2", "cv", "conf", "scores"),
    show_top: int = 10,
    output_dir: Optional[str] = None,
    plot: bool = False,
    split: bool = True,
    rank: str = "genus",
    use_names: bool = True,
    reord: str = "fpc",
) -> dict:
    """Mirror ``viz_annotation()`` in R: F1 / R2 / CV / confidence / scores plots."""
    results: dict = {}
    if data is None or data.empty:
        return results
    out_dir = Path(output_dir) if output_dir else None
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

    rank_key = "" if rank is None else str(rank).lower()
    skip_rank = rank_key in {"", "none", "exact", "taxid", "raw", "off", "false"}
    cache_file = None
    if out_dir and not skip_rank:
        cache_name = "taxid_genera_map.tsv" if rank_key in {"genus", "genera", "g"} else f"taxid_{rank_key}_map.tsv"
        cache_file = str(out_dir / cache_name)
    work = data.copy()
    if not skip_rank:
        work = remap_taxid_dataframe(work, rank=rank, cache_path=cache_file)
    rank_used = None if skip_rank else str(rank)

    selected = annotator_columns(work)
    work = _strip_annotator_prefixes(work)
    # After stripping, selected names lost prefixes
    annotators = []
    skip_meta = {"length", "sample", "true", "seq", "read_type"}
    for col in work.columns:
        if col in skip_meta or "confidence" in str(col).lower() or str(col).endswith("_conf"):
            continue
        annotators.append(col)

    taxid_xlab = "True taxon" if use_names else "True taxID"
    taxid_ylab = "Predicted taxon" if use_names else "Predicted taxID"
    if rank_used:
        taxid_xlab = f"{taxid_xlab} ({rank_used} or nearest rank)"
        taxid_ylab = f"{taxid_ylab} ({rank_used} or nearest rank)"
    italic = bool(use_names) and (rank_key in {"genus", "genera", "g", "species", "sp"} or skip_rank)

    name_map: Dict[str, str] = {}
    if use_names:
        ids = []
        for col in annotators + (["true"] if "true" in work.columns else []):
            ids.extend(work[col].dropna().astype(str).tolist())
        name_cache = str(out_dir / "taxid_names.tsv") if out_dir else None
        try:
            name_map = ensure_taxid_name_map(ids, cache_path=name_cache)
        except Exception:
            name_map = {}

    types = {str(t).lower() for t in type}
    from samovar.scores import (
        ENSEMBLE_NAME,
        has_usable_truth,
        is_platform_read_type,
        majority_vote,
        rows_with_usable_truth,
        tool_annotators,
    )

    truth_frame = rows_with_usable_truth(work) if has_usable_truth(work) else None
    platform_types = []
    if "read_type" in work.columns:
        platform_types = sorted(
            {
                str(x).strip().lower()
                for x in work["read_type"].fillna("")
                if is_platform_read_type(x)
            }
        )

    if types & {"scores", "score", "purity"}:
        from samovar.scores import (
            OPAL_DISPLAY,
            purity_by_taxon,
            save_scores_altair,
            save_scores_barplot,
            score_annotators,
        )

        scores_table = score_annotators(work, annotators)
        results["scores"] = scores_table
        if out_dir and not scores_table.empty:
            scores_table.to_csv(out_dir / "quality_scores.csv", index=False)
            caption_rank = f" ({rank_used})" if rank_used else ""
            save_scores_barplot(
                scores_table,
                out_dir / "scores.png",
                title=f"Annotation quality{caption_rank}",
            )
            save_scores_altair(scores_table, out_dir / "scores.html")
            save_scores_barplot(
                scores_table,
                out_dir / "opal_scores.png",
                title=f"OPAL-style profile metrics{caption_rank}",
                display=OPAL_DISPLAY,
                ymax=2.0,
            )
            save_scores_altair(
                scores_table,
                out_dir / "opal_scores.html",
                display=OPAL_DISPLAY,
            )
            if truth_frame is not None:
                taxon_frames = []
                tools = tool_annotators(annotators)
                seen = set()
                for name in scores_table["annotator"].astype(str):
                    if name in seen:
                        continue
                    seen.add(name)
                    if name == ENSEMBLE_NAME and name not in work.columns:
                        pred = majority_vote(work, tools).reindex(truth_frame.index)
                    elif name not in work.columns:
                        continue
                    else:
                        pred = work[name].reindex(truth_frame.index)
                    per = purity_by_taxon(truth_frame["true"], pred)
                    if per.empty:
                        continue
                    per.insert(0, "annotator", name)
                    taxon_frames.append(per)
                if taxon_frames:
                    pd.concat(taxon_frames, ignore_index=True).to_csv(
                        out_dir / "purity_by_taxon.csv", index=False
                    )
            try:
                from samovar.opal import maybe_run_opal

                plot_names = [
                    n for n in scores_table["annotator"].astype(str) if n in work.columns
                ]
                maybe_run_opal(
                    work,
                    plot_names,
                    out_dir,
                    rank=rank_used or "genus",
                )
            except Exception as exc:
                import logging

                logging.getLogger(__name__).warning("OPAL step skipped: %s", exc)

    if truth_frame is not None and types & {"f1"}:
        from samovar.scores import format_f1_caption, standard_classification_metrics

        gglist = {}
        true_levels_last: List[str] = []
        f1_work = truth_frame
        for name in annotators:
            pred = _collapse_other(f1_work[name], f1_work["true"])
            true = f1_work["true"].map(normalize_taxon_token)
            counts = pd.DataFrame({"true": true, "pred": pred}).value_counts().rename("Freq").reset_index()
            if reord == "fpc":
                fpc_core = fpc_taxon_order(counts, "pred", "true")
            else:
                fpc_core = sorted({t for t in list(counts["pred"]) + list(counts["true"]) if not is_special_taxon(t)})
            true_levels = axis_levels_from_fpc(true, fpc_core)
            pred_levels = axis_levels_from_fpc(pred, fpc_core)
            true_levels_last = true_levels
            f1 = _f1_from_counts(counts)
            std = standard_classification_metrics(true, pred)
            from samovar.opal import confusion_rates

            std.update(confusion_rates(true, pred))
            caption = format_f1_caption(f1, std, rank_used)
            matrix = _trim_matrix(_pivot_matrix(counts, "true", "pred", true_levels, pred_levels))
            gglist[name] = {"f1": f1, "metrics": std, "matrix": matrix}
            if out_dir:
                png = out_dir / f"F1_{name}.png"
                _save_heatmap_png(matrix, png, name, taxid_xlab, taxid_ylab, caption, name_map, italic)
                _save_altair(
                    _altair_heatmap(matrix, name, taxid_xlab, taxid_ylab),
                    out_dir / f"F1_{name}.html",
                    n_cells=int(matrix.size),
                )
                from samovar.stage_report import write_heatmap_mqc

                write_heatmap_mqc(
                    matrix,
                    out_dir / f"F1_{name}_mqc.json",
                    section_name=f"F1 heatmap — {name}",
                    description="",
                    xlab=taxid_xlab,
                    ylab=taxid_ylab,
                )
        results["F1"] = gglist
        if out_dir and len(platform_types) >= 2:
            from samovar.stage_report import write_heatmap_mqc

            for rt in platform_types:
                sub = f1_work[f1_work["read_type"].astype(str).str.strip().str.lower() == rt]
                if sub.empty or "read_type" not in f1_work.columns:
                    continue
                for name in annotators:
                    pred = _collapse_other(sub[name], sub["true"])
                    true = sub["true"].map(normalize_taxon_token)
                    counts = pd.DataFrame({"true": true, "pred": pred}).value_counts().rename("Freq").reset_index()
                    if counts.empty:
                        continue
                    if reord == "fpc":
                        fpc_core = fpc_taxon_order(counts, "pred", "true")
                    else:
                        fpc_core = sorted(
                            {t for t in list(counts["pred"]) + list(counts["true"]) if not is_special_taxon(t)}
                        )
                    matrix = _trim_matrix(
                        _pivot_matrix(
                            counts,
                            "true",
                            "pred",
                            axis_levels_from_fpc(true, fpc_core),
                            axis_levels_from_fpc(pred, fpc_core),
                        )
                    )
                    write_heatmap_mqc(
                        matrix,
                        out_dir / f"F1_{name}.{rt}_mqc.json",
                        section_name=f"F1 heatmap — {name}",
                        description="",
                        xlab=taxid_xlab,
                        ylab=taxid_ylab,
                    )

    if truth_frame is not None and types & {"r2"}:
        from samovar.scores import format_r2_caption

        gglist = {}
        r2_work = truth_frame
        for name in annotators:
            pred = _collapse_other(r2_work[name], r2_work["true"])
            true = r2_work["true"].map(normalize_taxon_token)
            table, r2 = _r2_table(true, pred)
            if np.isnan(r2):
                continue
            gglist[name] = {"r2": r2, "table": table}
            if out_dir:
                _save_scatter_png(table, out_dir / f"R2_{name}.png", name, format_r2_caption(r2, rank_used))
                _save_altair(_altair_scatter(table, name), out_dir / f"R2_{name}.html")
                from samovar.stage_report import write_scatter_mqc

                write_scatter_mqc(
                    table,
                    out_dir / f"R2_{name}_mqc.json",
                    section_name=f"Abundance R² — {name}",
                    description="",
                    xlab="Predicted taxon",
                    ylab="True taxon",
                )
        results["R2"] = gglist
        if out_dir and len(platform_types) >= 2 and "read_type" in r2_work.columns:
            from samovar.stage_report import write_scatter_mqc

            for rt in platform_types:
                sub = r2_work[r2_work["read_type"].astype(str).str.strip().str.lower() == rt]
                if sub.empty:
                    continue
                for name in annotators:
                    pred = _collapse_other(sub[name], sub["true"])
                    true = sub["true"].map(normalize_taxon_token)
                    table, r2 = _r2_table(true, pred)
                    if np.isnan(r2):
                        continue
                    write_scatter_mqc(
                        table,
                        out_dir / f"R2_{name}.{rt}_mqc.json",
                        section_name=f"Abundance R² — {name}",
                        description="",
                        xlab="Predicted taxon",
                        ylab="True taxon",
                    )

    if types & {"cv", "cross-validation"}:
        if len(annotators) >= 2:
            gglist = {}
            for i, name1 in enumerate(annotators):
                for name2 in annotators[:i]:
                    raw = pd.DataFrame(
                        {
                            name1: work[name1].map(normalize_taxon_token),
                            name2: work[name2].map(normalize_taxon_token),
                        }
                    )
                    if show_top != 0:
                        freq = pd.concat([raw[name1], raw[name2]]).value_counts()
                        keep = [t for t in freq.index if not is_special_taxon(t)][: max(1, show_top - 1)]
                        keep_set = set(keep) | {"0"}
                        raw[name1] = raw[name1].where(raw[name1].isin(keep_set), "other")
                        raw[name2] = raw[name2].where(raw[name2].isin(keep_set), "other")
                    counts = raw.value_counts().rename("Freq").reset_index()
                    counts = counts.rename(columns={name1: "a", name2: "b"})
                    if reord == "fpc":
                        fpc_core = fpc_taxon_order(counts, "b", "a")
                    else:
                        fpc_core = sorted({t for t in list(counts["a"]) + list(counts["b"]) if not is_special_taxon(t)})
                    lev1 = axis_levels_from_fpc(raw[name1], fpc_core)
                    lev2 = axis_levels_from_fpc(raw[name2], fpc_core)
                    matrix = _trim_matrix(
                        _pivot_matrix(
                            counts.rename(columns={"a": name1, "b": name2}),
                            name2,
                            name1,
                            lev2,
                            lev1,
                        )
                    )
                    key = f"{name1} vs {name2}"
                    gglist[key] = matrix
                    if out_dir:
                        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", key)
                        title = key if rank_used is None else f"{key} ({rank_used})"
                        _save_heatmap_png(
                            matrix,
                            out_dir / f"CV_{safe}.png",
                            title,
                            name2,
                            name1,
                            "",
                            name_map,
                            italic,
                        )
                        _save_altair(
                            _altair_heatmap(matrix, title, name2, name1),
                            out_dir / f"CV_{safe}.html",
                            n_cells=int(matrix.size),
                        )
                        from samovar.stage_report import write_heatmap_mqc

                        write_heatmap_mqc(
                            matrix,
                            out_dir / f"CV_{safe}_mqc.json",
                            section_name=f"Cross-validation — {title}",
                            description="",
                            xlab=name2,
                            ylab=name1,
                        )
            results["CV"] = gglist
            if out_dir and len(platform_types) >= 2:
                from samovar.stage_report import write_heatmap_mqc

                for rt in platform_types:
                    sub = work[work["read_type"].astype(str).str.strip().str.lower() == rt]
                    if sub.empty:
                        continue
                    for i, name1 in enumerate(annotators):
                        for name2 in annotators[:i]:
                            raw = pd.DataFrame(
                                {
                                    name1: sub[name1].map(normalize_taxon_token),
                                    name2: sub[name2].map(normalize_taxon_token),
                                }
                            )
                            counts = raw.value_counts().rename("Freq").reset_index()
                            if counts.empty:
                                continue
                            counts = counts.rename(columns={name1: "a", name2: "b"})
                            fpc_core = fpc_taxon_order(counts, "b", "a")
                            matrix = _trim_matrix(
                                _pivot_matrix(
                                    counts.rename(columns={"a": name1, "b": name2}),
                                    name2,
                                    name1,
                                    axis_levels_from_fpc(raw[name2], fpc_core),
                                    axis_levels_from_fpc(raw[name1], fpc_core),
                                )
                            )
                            safe = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{name1} vs {name2}")
                            write_heatmap_mqc(
                                matrix,
                                out_dir / f"CV_{safe}.{rt}_mqc.json",
                                section_name=f"Cross-validation — {name1} vs {name2}",
                                description="",
                                xlab=name2,
                                ylab=name1,
                            )

    conf_cols = [c for c in work.columns if str(c).endswith("_conf") or "_conf" in str(c)]
    if types & {"conf", "confidence"} and conf_cols:
        gglist = {}
        for conf_col in conf_cols:
            match = str(conf_col).replace("_conf", "")
            if match not in work.columns:
                continue
            tmp = pd.DataFrame({"x": pd.to_numeric(work[conf_col], errors="coerce"), "y": work[match].astype(str)})
            if show_top:
                top = tmp["y"].value_counts().head(show_top).index
                tmp = tmp[tmp["y"].isin(top)]
            gglist[match] = tmp
            if out_dir:
                _save_confidence_png(tmp, out_dir / f"confidence_{match}.png", f"{match} confidence")
        results["confidence"] = gglist

    return results


def compare_annotations(
    annotation_dir: str,
    output_dir: Optional[str] = None,
    csv_file: Optional[str] = None,
    show_top: int = 0,
    types: Sequence[str] = ("f1", "R2", "cv", "scores"),
    rank: str = "genus",
    split: bool = False,
) -> pd.DataFrame:
    """CLI analog of ``workflow/compare_annotations.R``."""
    require_viz_backend()
    annotation_path = Path(annotation_dir)
    if not annotation_path.is_dir():
        raise FileNotFoundError(f"annotation_dir not found: {annotation_path}")
    data = read_annotation_dir(annotation_dir)
    if data.empty:
        raise FileNotFoundError(
            f"no per-sample annotation tables in {annotation_path} "
            "(expected *.csv other than combined_annotation_table*)"
        )
    tax_cols = [c for c in data.columns if str(c).startswith("taxID_")]
    annotators = []
    for col in tax_cols:
        parts = str(col).split("_")
        if len(parts) >= 2:
            annotators.append(parts[1])
    type_list = list(types)
    if len(set(annotators)) < 2:
        type_list = [t for t in type_list if t.lower() not in {"cv", "cross-validation"}]
    viz_annotation(
        data,
        type=type_list,
        show_top=show_top,
        output_dir=output_dir,
        plot=False,
        split=split,
        rank=rank,
    )
    if csv_file:
        Path(csv_file).parent.mkdir(parents=True, exist_ok=True)
        data.to_csv(csv_file, index=False)
    return data
