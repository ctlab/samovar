"""Score regenerated abundance tables against the observed annotation.

The default builtin ``shannon_ks`` compares Shannon alpha-diversity
distributions (one value per annotator×sample) with a two-sample Kolmogorov–
Smirnov test. Lower KS statistic means the generated community better matches
the observed diversity; that candidate is selected for ISS.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

UNCLASSIFIED = frozenset({"0", "unclassified", "unclassified_root", "nan", "none", ""})
BUILTIN_TABLE_SCORERS = frozenset(
    {
        "shannon_ks",
        "ks",
        "shannon",
        "alpha_ks",
        "sparsedossa2_cv",
        "sparsedossa2-cv",
        "sparsedossa2",
        "sd2_cv",
        "sd2-cv",
        "fitcv",
    }
)

PathLike = Union[str, Path]


def shannon_entropy(counts: Iterable[float]) -> float:
    vec = np.asarray(list(counts), dtype=float)
    vec = vec[np.isfinite(vec) & (vec > 0)]
    total = float(vec.sum())
    if total <= 0:
        return 0.0
    p = vec / total
    return float(-(p * np.log(p)).sum())


def _drop_unclassified_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    if matrix.empty:
        return matrix
    keep = ~matrix.index.astype(str).str.lower().isin(UNCLASSIFIED)
    return matrix.loc[keep]


def shannon_vector_from_annotation(data: pd.DataFrame) -> np.ndarray:
    """Shannon diversity per (annotator, sample) from a long annotation table."""
    from samovar.regenerate import _count_matrix, _taxid_columns

    if data is None or data.empty:
        return np.array([], dtype=float)
    frame = data
    if "sample" not in frame.columns:
        frame = frame.copy()
        frame["sample"] = "1"
    values: List[float] = []
    for col in _taxid_columns(frame):
        matrix = _drop_unclassified_rows(_count_matrix(frame, col, "sample"))
        for sample in matrix.columns:
            values.append(shannon_entropy(matrix[sample].to_numpy()))
    return np.asarray(values, dtype=float)


def shannon_vector_from_table(table: pd.DataFrame) -> np.ndarray:
    values: List[float] = []
    if table is None or table.empty or "taxid" not in table.columns:
        return np.array([], dtype=float)
    work = table.copy()
    work["taxid"] = work["taxid"].astype(str)
    n_cols = [c for c in work.columns if str(c).startswith("N_")]
    work = work.loc[~work["taxid"].str.lower().isin(UNCLASSIFIED)]
    for col in n_cols:
        values.append(shannon_entropy(pd.to_numeric(work[col], errors="coerce").fillna(0)))
    return np.asarray(values, dtype=float)


def shannon_vector_from_tables(tables: Dict[str, pd.DataFrame]) -> np.ndarray:
    """Shannon diversity per (annotator table, N_* sample column)."""
    chunks = [shannon_vector_from_table(table) for table in (tables or {}).values()]
    nonempty = [c for c in chunks if c.size]
    if not nonempty:
        return np.array([], dtype=float)
    return np.concatenate(nonempty)


def shannon_vector_from_annotation_annotator(data: pd.DataFrame, annotator: str) -> np.ndarray:
    from samovar.regenerate import _annotator_name, _count_matrix, _taxid_columns

    if data is None or data.empty:
        return np.array([], dtype=float)
    frame = data
    if "sample" not in frame.columns:
        frame = frame.copy()
        frame["sample"] = "1"
    values: List[float] = []
    for col in _taxid_columns(frame):
        if _annotator_name(col) != annotator:
            continue
        matrix = _drop_unclassified_rows(_count_matrix(frame, col, "sample"))
        for sample in matrix.columns:
            values.append(shannon_entropy(matrix[sample].to_numpy()))
    return np.asarray(values, dtype=float)


def ks_shannon(
    observed: np.ndarray,
    generated: np.ndarray,
) -> Dict[str, Any]:
    """Two-sample KS test; ``rank_value`` is the statistic (lower is better)."""
    obs = np.asarray(observed, dtype=float)
    gen = np.asarray(generated, dtype=float)
    obs = obs[np.isfinite(obs)]
    gen = gen[np.isfinite(gen)]
    n_obs, n_gen = int(obs.size), int(gen.size)
    if n_obs == 0 or n_gen == 0:
        return _candidate_row(
            scorer="shannon_ks",
            rank_value=1.0,
            metric=1.0,
            metric_name="ks_statistic",
            pvalue=0.0,
            n_observed=n_obs,
            n_generated=n_gen,
            ks_statistic=1.0,
        )
    from scipy.stats import ks_2samp

    result = ks_2samp(obs, gen, alternative="two-sided", mode="auto")
    stat = float(result.statistic)
    pvalue = float(result.pvalue)
    if not math.isfinite(stat):
        stat = 1.0
    if not math.isfinite(pvalue):
        pvalue = 0.0
    return _candidate_row(
        scorer="shannon_ks",
        rank_value=stat,
        metric=stat,
        metric_name="ks_statistic",
        pvalue=pvalue,
        n_observed=n_obs,
        n_generated=n_gen,
        ks_statistic=stat,
        observed_shannon=[float(x) for x in obs],
        generated_shannon=[float(x) for x in gen],
    )


def _candidate_row(
    *,
    scorer: str,
    rank_value: float,
    metric: Optional[float] = None,
    metric_name: str = "rank_value",
    ok: bool = True,
    mode: Optional[str] = None,
    annotator: Optional[str] = None,
    pvalue: float = 0.0,
    n_observed: int = 0,
    n_generated: int = 0,
    error: Optional[str] = None,
    **extra: Any,
) -> Dict[str, Any]:
    """Shared table-scorer row: same keys for shannon_ks and sparsedossa2_cv."""
    row: Dict[str, Any] = {
        "scorer": scorer,
        "ok": bool(ok),
        "rank_value": float(rank_value) if math.isfinite(float(rank_value)) else 1.0e9,
        "metric": metric if metric is None or math.isfinite(float(metric)) else None,
        "metric_name": metric_name,
        "pvalue": float(pvalue) if pvalue is not None and math.isfinite(float(pvalue)) else 0.0,
        "n_observed": int(n_observed),
        "n_generated": int(n_generated),
    }
    if mode is not None:
        row["mode"] = mode
    if annotator is not None:
        row["annotator"] = annotator
    if error:
        row["error"] = str(error)[-2000:]
    if extra.get("ks_statistic") is not None:
        row["ks_statistic"] = extra["ks_statistic"]
    if extra.get("cv_goodness_of_fit") is not None:
        row["cv_goodness_of_fit"] = extra["cv_goodness_of_fit"]
    for key in ("observed_shannon", "generated_shannon", "details"):
        if key in extra and extra[key] is not None:
            row[key] = extra[key]
    return row


def canonicalize_table_scorer(name: Optional[str]) -> str:
    key = str(name or "shannon_ks").strip().lower().replace("-", "_")
    if key in {"sparsedossa2_cv", "sparsedossa2", "sd2_cv", "fitcv"}:
        return "sparsedossa2_cv"
    if key in {"shannon_ks", "ks", "shannon", "alpha_ks"}:
        return "shannon_ks"
    if not key:
        return "shannon_ks"
    raise ValueError(
        f"Unknown table scorer {name!r}. Built-in: shannon_ks, sparsedossa2_cv."
    )


def _tables_to_matrices(tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    matrices: Dict[str, pd.DataFrame] = {}
    for name, table in (tables or {}).items():
        if table is None or table.empty or "taxid" not in table.columns:
            continue
        work = table.copy()
        work["taxid"] = work["taxid"].astype(str)
        n_cols = [c for c in work.columns if str(c).startswith("N_")]
        if not n_cols:
            continue
        mat = work.set_index("taxid")[n_cols]
        mat.columns = [str(c)[2:] if str(c).startswith("N_") else str(c) for c in mat.columns]
        matrices[str(name)] = mat.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    return matrices


def score_sparsedossa2_cv(
    annotation: pd.DataFrame,
    tables: Dict[str, pd.DataFrame],
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """fitCV_SparseDOSSA2 on generated feature×sample tables. Higher GOF is better."""
    _ = annotation
    from samovar.sparsedossa2 import fitcv_score_matrix

    details: List[Dict[str, Any]] = []
    scores: List[float] = []
    for name, matrix in _tables_to_matrices(tables).items():
        payload = dict(fitcv_score_matrix(matrix, config=config or {}))
        payload["annotator"] = name
        details.append(payload)
        val = payload.get("cv_goodness_of_fit")
        try:
            num = float(val)
        except (TypeError, ValueError):
            num = float("nan")
        if math.isfinite(num):
            scores.append(num)
    if scores:
        mean = float(sum(scores) / len(scores))
        rank = -mean
    else:
        mean = float("nan")
        rank = 1.0e9
    n_gen = int(sum(m.shape[1] for m in _tables_to_matrices(tables).values()))
    n_obs = int(sum(m.shape[0] for m in _tables_to_matrices(tables).values()))
    return _candidate_row(
        scorer="sparsedossa2_cv",
        rank_value=rank,
        metric=mean if math.isfinite(mean) else None,
        metric_name="cv_goodness_of_fit",
        n_observed=n_obs,
        n_generated=n_gen,
        cv_goodness_of_fit=mean if math.isfinite(mean) else None,
        details=details,
    )


def _cv_cell_matrix(
    left: pd.DataFrame,
    right: pd.DataFrame,
    same: bool,
) -> pd.DataFrame:
    from samovar.sparsedossa2 import concat_aligned_matrices

    if same:
        return left
    return concat_aligned_matrices(left, right)


def score_sparsedossa2_cv_grid(
    tables_by_mode: Dict[str, Dict[str, pd.DataFrame]],
    annotator: str,
    modes: Sequence[str],
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """fitCV grid for one annotation table: methods × methods.

    Diagonal cells score that method's table. Off-diagonal cells score the
    sample-concatenated pair of two methods (same annotator only).
    """
    from samovar.sparsedossa2 import fitcv_score_jobs

    matrices: Dict[str, pd.DataFrame] = {}
    for mode in modes:
        tables = tables_by_mode.get(mode) or {}
        mats = _tables_to_matrices({annotator: tables[annotator]} if annotator in tables else {})
        if annotator in mats:
            matrices[mode] = mats[annotator]
    jobs: List[tuple] = []
    keys: List[tuple] = []
    for mode_i in modes:
        if mode_i not in matrices:
            continue
        for mode_j in modes:
            if mode_j not in matrices:
                continue
            key = f"{mode_i}__vs__{mode_j}"
            cell = _cv_cell_matrix(matrices[mode_i], matrices[mode_j], mode_i == mode_j)
            jobs.append((key, cell))
            keys.append((mode_i, mode_j, key))
    payloads = fitcv_score_jobs(jobs, config=config) if jobs else {}
    grid: Dict[str, Dict[str, Any]] = {m: {} for m in modes}
    diag_rows: List[Dict[str, Any]] = []
    for mode_i, mode_j, key in keys:
        payload = dict(payloads.get(key) or {})
        val = payload.get("cv_goodness_of_fit")
        try:
            num = float(val)
        except (TypeError, ValueError):
            num = float("nan")
        grid[mode_i][mode_j] = num
        if mode_i == mode_j:
            rank = -num if math.isfinite(num) else 1.0e9
            mat = matrices.get(mode_i)
            diag_rows.append(
                _candidate_row(
                    scorer="sparsedossa2_cv",
                    rank_value=rank,
                    metric=num if math.isfinite(num) else None,
                    metric_name="cv_goodness_of_fit",
                    ok=math.isfinite(num),
                    mode=mode_i,
                    annotator=annotator,
                    n_observed=int(mat.shape[0]) if mat is not None else 0,
                    n_generated=int(mat.shape[1]) if mat is not None else 0,
                    cv_goodness_of_fit=num if math.isfinite(num) else None,
                    error=payload.get("error"),
                    details=payload,
                )
            )
    winner = pick_best_table_method(diag_rows)
    if not winner:
        for mode in modes:
            if mode in matrices:
                winner = mode
                break
    return {
        "scorer": "sparsedossa2_cv",
        "annotator": annotator,
        "winner": winner,
        "score_matrix": grid,
        "cv_matrix": grid,
        "candidates": diag_rows,
        "metric_name": "cv_goodness_of_fit",
    }


def score_shannon_ks_annotator(
    annotation: pd.DataFrame,
    tables_by_mode: Dict[str, Dict[str, pd.DataFrame]],
    annotator: str,
    modes: Sequence[str],
) -> Dict[str, Any]:
    observed = shannon_vector_from_annotation_annotator(annotation, annotator)
    generated_by_mode: Dict[str, np.ndarray] = {}
    rows: List[Dict[str, Any]] = []
    for mode in modes:
        table = (tables_by_mode.get(mode) or {}).get(annotator)
        generated = shannon_vector_from_table(table) if table is not None else np.array([])
        generated_by_mode[mode] = generated
        score = ks_shannon(observed, generated)
        score["mode"] = mode
        score["annotator"] = annotator
        score["ok"] = bool(table is not None and not getattr(table, "empty", True))
        rows.append(score)
    grid: Dict[str, Dict[str, Any]] = {m: {} for m in modes}
    for mode_i in modes:
        for mode_j in modes:
            if mode_i == mode_j:
                grid[mode_i][mode_j] = next(
                    (r.get("ks_statistic") for r in rows if r.get("mode") == mode_i),
                    None,
                )
            else:
                pair = ks_shannon(generated_by_mode.get(mode_i, np.array([])), generated_by_mode.get(mode_j, np.array([])))
                grid[mode_i][mode_j] = pair.get("ks_statistic")
    return {
        "scorer": "shannon_ks",
        "annotator": annotator,
        "winner": pick_best_table_method(rows),
        "score_matrix": grid,
        "cv_matrix": grid,
        "candidates": rows,
        "metric_name": "ks_statistic",
    }


def rank_methods_per_annotator(
    annotation: pd.DataFrame,
    tables_by_mode: Dict[str, Dict[str, pd.DataFrame]],
    scorer: Optional[str] = None,
    config: Optional[Dict[str, Any]] = None,
    modes: Optional[Sequence[str]] = None,
) -> Dict[str, Any]:
    """Score regeneration methods independently for each annotation table."""
    kind = canonicalize_table_scorer(scorer)
    mode_list = list(modes or tables_by_mode.keys())
    annotators: List[str] = []
    seen = set()
    for tables in tables_by_mode.values():
        for name in tables or {}:
            if name not in seen:
                seen.add(name)
                annotators.append(name)
    by_annotator: Dict[str, Any] = {}
    winner_by_annotator: Dict[str, str] = {}
    mixed: Dict[str, pd.DataFrame] = {}
    for annotator in annotators:
        if kind == "sparsedossa2_cv":
            block = score_sparsedossa2_cv_grid(
                tables_by_mode, annotator, mode_list, config=config
            )
        else:
            block = score_shannon_ks_annotator(
                annotation, tables_by_mode, annotator, mode_list
            )
        by_annotator[annotator] = block
        winner = block.get("winner")
        if winner:
            winner_by_annotator[annotator] = str(winner)
            table = (tables_by_mode.get(winner) or {}).get(annotator)
            if table is not None and not table.empty:
                mixed[annotator] = table
        if annotator not in mixed:
            for mode in mode_list:
                table = (tables_by_mode.get(mode) or {}).get(annotator)
                if table is not None and not table.empty:
                    mixed[annotator] = table
                    winner_by_annotator.setdefault(annotator, str(mode))
                    break
    winners = list(winner_by_annotator.values())
    overall = winners[0] if winners and all(w == winners[0] for w in winners) else "mixed"
    flat = []
    for block in by_annotator.values():
        flat.extend(block.get("candidates") or [])
    metric_name = "rank_value"
    for block in by_annotator.values():
        if block.get("metric_name"):
            metric_name = str(block["metric_name"])
            break
    return {
        "scorer": kind,
        "winner": overall,
        "winner_by_annotator": winner_by_annotator,
        "by_annotator": by_annotator,
        "candidates": flat,
        "tables": mixed,
        "metric_name": metric_name,
    }


def score_generated_tables(
    annotation: pd.DataFrame,
    tables: Dict[str, pd.DataFrame],
    scorer: Optional[str] = None,
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    kind = canonicalize_table_scorer(scorer)
    if kind == "shannon_ks":
        return ks_shannon(
            shannon_vector_from_annotation(annotation),
            shannon_vector_from_tables(tables),
        )
    if kind == "sparsedossa2_cv":
        return score_sparsedossa2_cv(annotation, tables, config=config)
    raise ValueError(f"Unsupported table scorer {kind!r}")


def pick_best_table_method(rows: Sequence[Dict[str, Any]]) -> Optional[str]:
    """Smallest KS statistic, then largest p-value, then original order."""
    ranked = []
    for index, row in enumerate(rows):
        if not row.get("ok", True):
            continue
        ranked.append(
            (
                float(row.get("rank_value", 1.0)),
                -float(row.get("pvalue", 0.0)),
                index,
                str(row.get("mode") or ""),
            )
        )
    if not ranked:
        return None
    ranked.sort()
    return ranked[0][3]


def write_table_selection(
    dest: PathLike,
    payload: Dict[str, Any],
) -> Path:
    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)

    def _sanitize(obj: Any) -> Any:
        if isinstance(obj, dict):
            return {str(k): _sanitize(v) for k, v in obj.items()}
        if isinstance(obj, (list, tuple)):
            return [_sanitize(v) for v in obj]
        if isinstance(obj, float) and not math.isfinite(obj):
            return None
        if isinstance(obj, (np.floating, np.integer)):
            val = obj.item()
            if isinstance(val, float) and not math.isfinite(val):
                return None
            return val
        if isinstance(obj, np.ndarray):
            return _sanitize(obj.tolist())
        return obj

    path.write_text(json.dumps(_sanitize(payload), indent=2) + "\n", encoding="utf-8")
    return path


CANDIDATE_KEYS = (
    "scorer",
    "ok",
    "rank_value",
    "metric",
    "metric_name",
    "pvalue",
    "n_observed",
    "n_generated",
)


def _viz_backend_wanted(name: str, config: Optional[Dict[str, Any]] = None) -> bool:
    selected = (config or {}).get("scoring_tools")
    if selected is None:
        return True
    tokens = {str(x).strip().lower().replace("-", "_") for x in selected if str(x).strip()}
    if not tokens:
        return True
    viz = {"altair", "cnsplots", "matplotlib"}
    if tokens & viz:
        return name.lower().replace("-", "_") in tokens
    return True


def _grid_to_frame(grid: Dict[str, Dict[str, Any]], modes: Sequence[str]) -> pd.DataFrame:
    rows = []
    for mode_i in modes:
        row = []
        for mode_j in modes:
            val = (grid.get(mode_i) or {}).get(mode_j)
            try:
                num = float(val)
            except (TypeError, ValueError):
                num = float("nan")
            row.append(num if math.isfinite(num) else np.nan)
        rows.append(row)
    return pd.DataFrame(rows, index=list(modes), columns=list(modes))


def table_score_plot_dirs(abundance_dir: PathLike) -> List[Path]:
    out = Path(abundance_dir)
    dirs = [out / "table_score_plots"]
    if out.name == ".regenerated_abundance":
        parent = out.parent
        run = parent.parent if parent.name == "regenerated" else parent
        dirs.append(run / "regenerated_annotations_plots")
    return dirs


def _save_table_score_heatmap_png(matrix: pd.DataFrame, path: Path, title: str, xlabel: str, ylabel: str) -> None:
    try:
        from samovar.viz_annotation import _setup_cns, _use_agg_backend
    except Exception:
        return
    try:
        _use_agg_backend()
        cns = _setup_cns()
        import matplotlib.pyplot as plt
        from matplotlib.colors import LinearSegmentedColormap
    except Exception:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    values = matrix.to_numpy(dtype=float)
    n = max(matrix.shape[0], matrix.shape[1], 1)
    fig, ax = plt.subplots(figsize=(max(4.5, 0.7 * n + 2), max(4.0, 0.7 * n + 1.8)), facecolor="white")
    ax.set_facecolor("white")
    cmap = LinearSegmentedColormap.from_list("samovar_table_score", ["#FFFFFF", "#78C679", "#004529"])
    masked = np.ma.masked_invalid(values)
    im = ax.imshow(masked, cmap=cmap, aspect="equal", origin="upper")
    ax.set_xticks(range(len(matrix.columns)))
    ax.set_yticks(range(len(matrix.index)))
    ax.set_xticklabels([str(c) for c in matrix.columns], rotation=35, ha="right")
    ax.set_yticklabels([str(c) for c in matrix.index])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    fig.tight_layout()
    if cns is not None and hasattr(cns, "savefig"):
        cns.savefig(str(path))
    else:
        fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close("all")


def write_table_score_plots(
    payload: Dict[str, Any],
    dest: PathLike,
    config: Optional[Dict[str, Any]] = None,
) -> List[Path]:
    """MultiQC JSON plus optional cnsplots PNG / Altair HTML for table scoring."""
    from samovar.scores import SCORE_DISPLAY, save_scores_altair, save_scores_barplot
    from samovar.stage_report import write_bargraph_mqc, write_heatmap_mqc, write_table_mqc

    out = Path(dest)
    out.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []
    by_annotator = payload.get("by_annotator") or {}
    if not by_annotator and payload.get("candidates"):
        by_annotator = {"all": {"candidates": payload["candidates"], "score_matrix": {}}}
    want_png = _viz_backend_wanted("cnsplots", config) or _viz_backend_wanted("matplotlib", config)
    want_html = _viz_backend_wanted("altair", config)
    try:
        import altair  # noqa: F401
    except ImportError:
        want_html = False

    all_rows: List[Dict[str, Any]] = []
    for annotator, block in by_annotator.items():
        rows = list(block.get("candidates") or [])
        for row in rows:
            labeled = dict(row)
            labeled["label"] = f"{annotator} / {row.get('mode') or 'method'}"
            all_rows.append(labeled)
        metric_name = str(block.get("metric_name") or payload.get("metric_name") or "rank_value")
        metric_label = SCORE_DISPLAY.get(metric_name, metric_name.replace("_", " "))
        grid = block.get("score_matrix") or block.get("cv_matrix") or {}
        modes = list(grid.keys()) if grid else [str(r.get("mode")) for r in rows if r.get("mode")]
        if grid and modes:
            matrix = _grid_to_frame(grid, modes)
            slug = re.sub(r"[^A-Za-z0-9._-]+", "_", str(annotator)).strip("_") or "annotator"
            mqc = out / f"TableScore_{slug}_mqc.json"
            finite = [v for v in matrix.to_numpy(dtype=float).ravel() if math.isfinite(v)]
            min_value = 0.0 if metric_name == "ks_statistic" else None
            if min_value is None and finite and min(finite) >= 0:
                min_value = 0.0
            write_heatmap_mqc(
                matrix,
                mqc,
                section_name=f"Table score — {annotator}",
                description=f"{metric_label} between table_reads_generator methods (diagonal vs observed).",
                xlab="Method",
                ylab="Method",
                min_value=min_value,
            )
            written.append(mqc)
            if want_png:
                png = out / f"TableScore_{slug}.png"
                try:
                    _save_table_score_heatmap_png(
                        matrix,
                        png,
                        f"Table score — {annotator}",
                        "Method",
                        "Method",
                    )
                    if png.is_file():
                        written.append(png)
                except Exception:
                    pass
            if want_html:
                try:
                    from samovar.viz_annotation import _altair_heatmap, _save_altair

                    html = out / f"TableScore_{slug}.html"
                    _save_altair(
                        _altair_heatmap(matrix, f"Table score — {annotator}", "Method", "Method"),
                        html,
                        n_cells=int(matrix.size),
                    )
                    if html.is_file() or html.with_suffix(".json").is_file():
                        written.append(html if html.is_file() else html.with_suffix(".json"))
                except Exception:
                    pass
        if rows:
            series: Dict[str, Dict[str, float]] = {metric_label: {}}
            for row in rows:
                mode = str(row.get("mode") or "method")
                val = row.get(metric_name)
                if val is None:
                    val = row.get("metric")
                if val is None or not math.isfinite(float(val)):
                    continue
                series[metric_label][mode] = float(val)
            if series[metric_label]:
                bars = out / f"TableScore_{re.sub(r'[^A-Za-z0-9._-]+', '_', str(annotator)).strip('_')}_bars_mqc.json"
                write_bargraph_mqc(
                    series,
                    bars,
                    section_name=f"Table score bars — {annotator}",
                    description=metric_label,
                    xlab=metric_label,
                )
                written.append(bars)
            frame = pd.DataFrame(rows)
            if not frame.empty and "mode" in frame.columns:
                plot_frame = frame.rename(columns={"mode": "annotator"})
                display = {metric_name: metric_label}
                if want_png:
                    try:
                        save_scores_barplot(
                            plot_frame,
                            out / f"TableScore_{re.sub(r'[^A-Za-z0-9._-]+', '_', str(annotator)).strip('_')}_bars.png",
                            title=f"Table score — {annotator}",
                            display=display,
                        )
                    except Exception:
                        pass
                if want_html:
                    try:
                        save_scores_altair(
                            plot_frame,
                            out / f"TableScore_{re.sub(r'[^A-Za-z0-9._-]+', '_', str(annotator)).strip('_')}_bars.html",
                            display=display,
                        )
                    except Exception:
                        pass

    if all_rows:
        table_path = out / "TableScore_quality_scores_mqc.json"
        write_table_mqc(
            all_rows,
            table_path,
            section_name="Table generator scores",
            description="Per-annotator ranking of table_reads_generator methods.",
            col1_header="Annotator / method",
            id_field="label",
            numeric_fields=(
                "metric",
                "rank_value",
                "pvalue",
                "ks_statistic",
                "cv_goodness_of_fit",
                "n_observed",
                "n_generated",
            ),
        )
        written.append(table_path)
    return written
