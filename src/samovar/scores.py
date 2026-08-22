"""Quality scores for taxonomic annotation vs a known true taxid.

Accuracy purity asks: among reads that truly come from one taxid, how much of
that mass did the tool pile onto a single predicted label (the majority call)?
F1 purity keeps that recall and penalises other true taxa that were also
assigned to the same majority label.
"""

from __future__ import annotations

from collections import Counter
import re
from typing import Dict, Iterable, List, Optional, Sequence

import numpy as np
import pandas as pd


def _taxon_helpers():
    from samovar.viz_annotation import is_special_taxon, normalize_taxon_token

    return is_special_taxon, normalize_taxon_token

COMBINATION_ANNOTATORS = frozenset({"samovar", "consensus"})

SCORE_COLUMNS = (
    "annotator",
    "n_reads",
    "n_taxa",
    "accuracy",
    "precision_macro",
    "recall_macro",
    "f1_micro",
    "f1_macro",
    "f1_weighted",
    "f1",
    "r2",
    "accuracy_purity",
    "f1_purity",
)


def majority_vote(frame: pd.DataFrame, annotators: Sequence[str]) -> pd.Series:
    """Per-read consensus: mode of non-unclassified votes, else unclassified."""
    is_special_taxon, normalize_taxon_token = _taxon_helpers()
    if not annotators:
        return pd.Series(dtype=str)
    cols = [c for c in annotators if c in frame.columns]
    if not cols:
        return pd.Series(dtype=str)

    def _row_vote(row) -> str:
        tokens = [normalize_taxon_token(row[c]) for c in cols]
        classified = [t for t in tokens if not is_special_taxon(t)]
        pool = classified or tokens
        counts = Counter(pool)
        top = counts.most_common()
        if not top:
            return "0"
        best_n = top[0][1]
        tied = [lab for lab, n in top if n == best_n]
        if len(tied) == 1:
            return tied[0]
        for token in tokens:
            if token in tied:
                return token
        return tied[0]

    return frame[cols].apply(_row_vote, axis=1)


def _confusion(true: Iterable, pred: Iterable) -> pd.DataFrame:
    _, normalize_taxon_token = _taxon_helpers()
    tmp = pd.DataFrame(
        {
            "true": pd.Series(true, dtype=object).map(normalize_taxon_token),
            "pred": pd.Series(pred, dtype=object).map(normalize_taxon_token),
        }
    )
    if tmp.empty:
        return pd.DataFrame(index=pd.Index([], name="true"), columns=pd.Index([], name="pred"))
    return pd.crosstab(tmp["true"], tmp["pred"]).astype(float)


def _f1_accuracy(true: pd.Series, pred: pd.Series) -> float:
    """Current heatmap F1: P = R = TP/N, so F1 equals per-read accuracy."""
    _, normalize_taxon_token = _taxon_helpers()
    t = true.map(normalize_taxon_token)
    p = pred.map(normalize_taxon_token)
    n = len(t)
    if n == 0:
        return 0.0
    tp = float((t == p).sum())
    return float(tp / n)


def standard_classification_metrics(true: Iterable, pred: Iterable) -> Dict[str, float]:
    """sklearn-style single-label scores (macro / weighted / micro / accuracy)."""
    from sklearn.metrics import (
        accuracy_score,
        f1_score,
        precision_score,
        recall_score,
    )

    _, normalize_taxon_token = _taxon_helpers()
    y = pd.Series(list(true) if not isinstance(true, pd.Series) else true).map(
        normalize_taxon_token
    ).astype(str)
    p = pd.Series(list(pred) if not isinstance(pred, pd.Series) else pred).map(
        normalize_taxon_token
    ).astype(str)
    if y.empty:
        return {
            "accuracy": 0.0,
            "precision_macro": 0.0,
            "recall_macro": 0.0,
            "f1_micro": 0.0,
            "f1_macro": 0.0,
            "f1_weighted": 0.0,
        }
    labels = sorted(set(y) | set(p))
    kwargs = {"labels": labels, "zero_division": 0}
    return {
        "accuracy": float(accuracy_score(y, p)),
        "precision_macro": float(precision_score(y, p, average="macro", **kwargs)),
        "recall_macro": float(recall_score(y, p, average="macro", **kwargs)),
        "f1_micro": float(f1_score(y, p, average="micro", **kwargs)),
        "f1_macro": float(f1_score(y, p, average="macro", **kwargs)),
        "f1_weighted": float(f1_score(y, p, average="weighted", **kwargs)),
    }


def _r2_abundance(true: pd.Series, pred: pd.Series) -> float:
    _, normalize_taxon_token = _taxon_helpers()
    t = true.map(normalize_taxon_token)
    p = pred.map(normalize_taxon_token)
    true_n = t.value_counts()
    pred_n = p.value_counts()
    taxa = sorted(set(true_n.index) | set(pred_n.index))
    if not taxa:
        return float("nan")
    y = np.array([float(true_n.get(x, 0)) for x in taxa])
    x = np.array([float(pred_n.get(x, 0)) for x in taxa])
    ss_tot = float(((y - y.mean()) ** 2).sum())
    ss_res = float(((x - y) ** 2).sum())
    if ss_tot == 0:
        return 0.0
    return float(1.0 - ss_res / ss_tot)


def purity_by_taxon(true: Iterable, pred: Iterable) -> pd.DataFrame:
    """Per-true-taxid majority recall / precision / F1.

    For true taxid ``t`` let ``m(t)`` be the most abundant *classified*
    prediction among reads whose true label is only ``t``. Then

    * recall = n(true=t, pred=m(t)) / n(true=t)
    * precision = n(true=t, pred=m(t)) / n(pred=m(t))
    * F1 = harmonic mean of precision and recall
    """
    is_special_taxon, _ = _taxon_helpers()
    mat = _confusion(true, pred)
    if mat.empty:
        return pd.DataFrame(
            columns=[
                "true",
                "majority_pred",
                "n_true",
                "n_majority",
                "n_predicted_as_majority",
                "recall",
                "precision",
                "f1",
            ]
        )
    pred_totals = mat.sum(axis=0)
    rows = []
    for taxid in mat.index:
        if is_special_taxon(taxid):
            continue
        row = mat.loc[taxid]
        n_true = float(row.sum())
        if n_true <= 0:
            continue
        classified = {
            str(lab): float(val)
            for lab, val in row.items()
            if not is_special_taxon(lab)
        }
        if classified:
            majority = max(classified.items(), key=lambda kv: (kv[1], kv[0]))[0]
            n_maj = classified[majority]
        else:
            majority = "0"
            n_maj = 0.0
        n_pred = float(pred_totals.get(majority, 0.0)) if majority != "0" else 0.0
        recall = n_maj / n_true
        precision = n_maj / n_pred if n_pred > 0 else 0.0
        if precision + recall == 0:
            f1 = 0.0
        else:
            f1 = 2 * precision * recall / (precision + recall)
        rows.append(
            {
                "true": str(taxid),
                "majority_pred": str(majority),
                "n_true": n_true,
                "n_majority": n_maj,
                "n_predicted_as_majority": n_pred,
                "recall": recall,
                "precision": precision,
                "f1": f1,
            }
        )
    return pd.DataFrame(rows)


def accuracy_purity(true: Iterable, pred: Iterable) -> float:
    """Read-weighted majority recall (cluster purity of true taxa).

    ``sum_t n(t, m(t)) / sum_t n(t)`` where ``m(t)`` is the most abundant
    classified prediction among reads from true taxid ``t`` only.
    """
    table = purity_by_taxon(true, pred)
    if table.empty:
        return 0.0
    denom = float(table["n_true"].sum())
    if denom <= 0:
        return 0.0
    return float(table["n_majority"].sum() / denom)


def f1_purity(true: Iterable, pred: Iterable) -> float:
    """Support-weighted F1 of per-taxon majority mapping.

    Recall is accuracy purity of taxid ``t``. Precision divides that majority
    count by every read the tool assigned to ``m(t)``, including reads from
    other true taxa.
    """
    table = purity_by_taxon(true, pred)
    if table.empty:
        return 0.0
    weights = table["n_true"].to_numpy(dtype=float)
    denom = float(weights.sum())
    if denom <= 0:
        return 0.0
    return float(np.dot(table["f1"].to_numpy(dtype=float), weights) / denom)


def annotation_scores(true: Iterable, pred: Iterable) -> Dict[str, float]:
    true_s = pd.Series(list(true) if not isinstance(true, pd.Series) else true.reset_index(drop=True))
    pred_s = pd.Series(list(pred) if not isinstance(pred, pd.Series) else pred.reset_index(drop=True))
    n = len(true_s)
    table = purity_by_taxon(true_s, pred_s)
    stats = {
        "n_reads": float(n),
        "n_taxa": float(len(table)),
        "accuracy_purity": accuracy_purity(true_s, pred_s),
        "f1_purity": f1_purity(true_s, pred_s),
        "f1": _f1_accuracy(true_s, pred_s),
        "r2": _r2_abundance(true_s, pred_s),
    }
    stats.update(standard_classification_metrics(true_s, pred_s))
    return stats


def _is_combination(name: str) -> bool:
    token = re.sub(r"^taxid[_./]*", "", str(name).strip(), flags=re.I).lower()
    return token in COMBINATION_ANNOTATORS


def tool_annotators(names: Sequence[str]) -> List[str]:
    return [n for n in names if not _is_combination(n)]


def order_annotators(names: Sequence[str]) -> List[str]:
    tools = sorted(tool_annotators(names), key=lambda n: str(n).lower())
    combo = [n for n in names if _is_combination(n)]
    combo.sort(key=lambda n: (0 if str(n).lower() == "consensus" else 1, str(n).lower()))
    return tools + combo


def score_annotators(
    work: pd.DataFrame,
    annotators: Sequence[str],
    true_col: str = "true",
    include_consensus: bool = True,
) -> pd.DataFrame:
    """Score each annotator plus optional majority-vote consensus / SAMOVAR."""
    if true_col not in work.columns:
        return pd.DataFrame(columns=list(SCORE_COLUMNS))
    true = work[true_col]
    names = list(annotators)
    tools = tool_annotators(names)
    rows = []

    def _add(name: str, pred: pd.Series) -> None:
        from samovar.viz_annotation import _collapse_other

        collapsed = _collapse_other(pred, true)
        stats = annotation_scores(true, pred)
        # Heatmap F1 / R² / sklearn metrics use the same collapsed labels as the plots.
        stats["f1"] = _f1_accuracy(true, collapsed)
        stats["r2"] = _r2_abundance(true, collapsed)
        stats.update(standard_classification_metrics(true, collapsed))
        rows.append({"annotator": name, **stats})

    for name in names:
        if name not in work.columns:
            continue
        _add(name, work[name])

    have = {r["annotator"] for r in rows}
    if include_consensus and len(tools) >= 2 and "consensus" not in have:
        _add("consensus", majority_vote(work, tools))

    table = pd.DataFrame(rows)
    if table.empty:
        return pd.DataFrame(columns=list(SCORE_COLUMNS))
    table = table[list(SCORE_COLUMNS)]
    ordered = order_annotators(list(table["annotator"]))
    table["annotator"] = pd.Categorical(table["annotator"], categories=ordered, ordered=True)
    return table.sort_values("annotator").reset_index(drop=True)


SCORE_DISPLAY = {
    "accuracy": "Accuracy",
    "precision_macro": "P-macro",
    "recall_macro": "R-macro",
    "f1_macro": "F1-macro",
    "f1": "F1 (current)",
    "r2": "R² (current)",
    "accuracy_purity": "Accuracy purity",
    "f1_purity": "F1 purity",
}

BAR_COLORS = {
    "accuracy": "#66A61E",
    "precision_macro": "#E6AB02",
    "recall_macro": "#A6761D",
    "f1_macro": "#1F78B4",
    "f1": "#7570B3",
    "r2": "#E7298A",
    "accuracy_purity": "#1B9E77",
    "f1_purity": "#D95F02",
}


def scores_long(table: pd.DataFrame) -> pd.DataFrame:
    metrics = [c for c in SCORE_DISPLAY if c in table.columns]
    long = table.melt(
        id_vars=["annotator"],
        value_vars=metrics,
        var_name="metric",
        value_name="score",
    )
    long["label"] = long["metric"].map(SCORE_DISPLAY)
    return long


def save_scores_barplot(table: pd.DataFrame, path, title: Optional[str] = None) -> None:
    """Grouped barplot of quality metrics per annotator (PNG)."""
    from pathlib import Path

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    from samovar.viz_annotation import _setup_cns

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if table is None or table.empty:
        return

    cns = _setup_cns()
    metrics = [m for m in SCORE_DISPLAY if m in table.columns]
    annotators = [str(a) for a in table["annotator"].tolist()]
    n_ann = len(annotators)
    n_met = len(metrics)
    x = np.arange(n_ann, dtype=float)
    width = 0.8 / max(n_met, 1)

    if cns is not None and hasattr(cns, "figure"):
        cns.figure(max(220, 48 * n_ann + 90), 190)
        ax = plt.gca()
    else:
        fig_w = max(7.2, 1.15 * n_ann + 2.4)
        _, ax = plt.subplots(figsize=(fig_w, 5.2))

    for i, metric in enumerate(metrics):
        vals = pd.to_numeric(table[metric], errors="coerce").to_numpy(dtype=float)
        offset = (i - (n_met - 1) / 2) * width
        ax.bar(
            x + offset,
            vals,
            width,
            label=SCORE_DISPLAY[metric],
            color=BAR_COLORS.get(metric, "0.5"),
            edgecolor="white",
            linewidth=0.4,
        )

    ax.set_xticks(x)
    ax.set_xticklabels(annotators, rotation=25, ha="right")
    ax.set_ylabel("Score")
    ax.set_xlabel("Annotator")
    ax.set_title(title or "Annotation quality")
    finite = pd.to_numeric(table[metrics].stack(), errors="coerce")
    ymax = float(np.nanmax(finite.to_numpy())) if len(finite) else 1.0
    ymin = float(np.nanmin(finite.to_numpy())) if len(finite) else 0.0
    ax.set_ylim(min(0.0, ymin) - 0.05, max(1.0, ymax) + 0.08)
    ax.axhline(0, color="0.6", linewidth=0.6)
    ax.legend(frameon=False, ncol=min(n_met, 4), loc="upper center", bbox_to_anchor=(0.5, -0.28))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True, linestyle=":", linewidth=0.5, color="0.8")
    fig = ax.figure
    fig.tight_layout()
    fig.savefig(
        path,
        dpi=150,
        bbox_inches="tight",
        pad_inches=0.35,
        facecolor="white",
        edgecolor="none",
    )
    plt.close("all")


def format_f1_caption(f1: float, std: Dict[str, float], rank: Optional[str] = None) -> str:
    """Two-line caption for F1 heatmaps: current formula + standard metrics."""
    tag = f" ({rank})" if rank else ""
    acc = std.get("accuracy", f1)
    return (
        f"F1 (current){tag} = 2PR/(P+R) = {f1:.3f}   |   P = R = TP/N = accuracy = {acc:.3f}\n"
        f"P-macro = {std.get('precision_macro', float('nan')):.3f}   "
        f"R-macro = {std.get('recall_macro', float('nan')):.3f}   "
        f"F1-macro = {std.get('f1_macro', float('nan')):.3f}   "
        f"F1-weighted = {std.get('f1_weighted', float('nan')):.3f}"
    )


def format_r2_caption(r2: float, rank: Optional[str] = None) -> str:
    """Two-line caption for abundance R² scatter plots (current formula)."""
    tag = f" ({rank})" if rank else ""
    return (
        f"R² (current){tag} = 1 − Σ(n̂_t − n_t)² / Σ(n_t − n̄)² = {r2:.3f}\n"
        f"n_t = true read count of taxid t; n̂_t = predicted count (labels not in truth → other)"
    )


def save_scores_altair(table: pd.DataFrame, path) -> None:
    from pathlib import Path

    path = Path(path)
    if table is None or table.empty:
        return
    try:
        import altair as alt
    except ImportError:
        return
    long = scores_long(table)
    chart = (
        alt.Chart(long)
        .mark_bar()
        .encode(
            x=alt.X("annotator:N", title="Annotator", sort=list(table["annotator"].astype(str))),
            y=alt.Y("score:Q", title="Score"),
            color=alt.Color("label:N", title="Metric"),
            xOffset="label:N",
            tooltip=["annotator", "label", "score"],
        )
        .properties(title="Annotation quality")
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        chart.save(str(path))
    except Exception:
        try:
            path.with_suffix(".json").write_text(chart.to_json())
        except Exception:
            pass
