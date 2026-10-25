"""Optional CAMI OPAL integration and OPAL-style profiling metrics.

OPAL (https://github.com/CAMI-challenge/OPAL) scores *taxonomic profiles*
(presence and relative abundance of taxa), not per-read labels. SamovaR:

* always computes the same presence / abundance metrics in Python
* if ``opal.py`` is installed (``./install.sh OPAL``), also writes CAMI
  ``.profile`` files and runs OPAL into ``<plots>/opal/``
"""

from __future__ import annotations

import logging
import os
import subprocess
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

CAMI_RANK_ALIASES = {
    "g": "genus",
    "genera": "genus",
    "sp": "species",
    "s": "species",
    "none": "strain",
    "exact": "strain",
    "taxid": "strain",
    "raw": "strain",
    "off": "strain",
    "false": "strain",
}
OPAL_RANKS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
    "strain",
)


def _taxon_helpers():
    from samovar.viz_annotation import is_special_taxon, normalize_taxon_token  # noqa: F401

    return is_special_taxon, normalize_taxon_token


def normalize_series(values):
    from samovar.viz_annotation import normalize_series as _f

    return _f(values)


def special_mask(values):
    from samovar.viz_annotation import special_mask as _f

    return _f(values)


def cami_rank(rank: Optional[str]) -> str:
    key = "" if rank is None else str(rank).strip().lower()
    mapped = CAMI_RANK_ALIASES.get(key, key or "genus")
    return mapped if mapped in OPAL_RANKS else "genus"


def opal_rank_range(rank: Optional[str]) -> str:
    """OPAL ``-r lo,hi``: include finer ranks so phages without a genus still score."""
    requested = cami_rank(rank)
    order = list(OPAL_RANKS)
    if requested not in order:
        return f"{requested},{requested}"
    start = order.index(requested)
    end = order.index("species") if "species" in order else start
    lo, hi = (start, end) if start <= end else (end, start)
    if lo == hi:
        lo = max(0, lo - 1)
    return f"{order[lo]},{order[hi]}"


def opal_enabled() -> bool:
    flag = os.environ.get("SAMOVAR_OPAL", "").strip().lower()
    if flag in {"0", "false", "no", "off"}:
        return False
    try:
        from samovar.paths import load_config

        cfg = load_config()
        if str(cfg.get("opal") or "").strip().lower() in {"0", "false", "no"}:
            return False
    except Exception:
        pass
    return True


def opal_executable() -> Optional[str]:
    """Path to ``opal.py`` / ``opal`` if the optional CAMI tool is installed."""
    from samovar.paths import discover_opal

    found = discover_opal()
    if found:
        return found
    return None


def opal_command() -> Optional[List[str]]:
    """Argv prefix to invoke OPAL (script, or ``python opal.py`` if not executable)."""
    exe = opal_executable()
    if not exe:
        return None
    if exe.endswith(".py") or not os.access(exe, os.X_OK):
        from samovar.paths import python_path

        return [python_path(), "-m", "samovar.opal", exe]
    return [exe]


def write_cami_profile(
    counts: Dict[str, float],
    dest: Path,
    sample_id: str = "1",
    rank: str = "genus",
    total: Optional[float] = None,
) -> Path:
    """Write a single-sample CAMI bioboxes taxonomic profile (RFC 0.10.0)."""
    from samovar.cami_profile import write_cami_profile as _write

    return _write(counts, dest, sample_id=sample_id, rank=rank, total=total)


def series_to_counts(values: Iterable) -> Dict[str, float]:
    is_special_taxon, normalize_taxon_token = _taxon_helpers()
    counts: Dict[str, float] = {}
    for raw in values:
        token = normalize_taxon_token(raw)
        if is_special_taxon(token):
            continue
        counts[token] = counts.get(token, 0.0) + 1.0
    return counts


def confusion_rates(true: Iterable, pred: Iterable) -> Dict[str, float]:
    """Macro one-vs-rest rates, including true-negative and false-negative %."""
    is_special_taxon, normalize_taxon_token = _taxon_helpers()
    y = normalize_series(true)
    p = normalize_series(pred)
    classes = [t for t in pd.unique(y) if not is_special_taxon(t)]
    empty = {
        "tpr": 0.0,
        "tnr": 0.0,
        "fpr": 0.0,
        "fnr": 0.0,
        "tn_pct": 0.0,
        "fn_pct": 0.0,
        "fp_pct": 0.0,
        "tp_pct": 0.0,
    }
    if y.empty or not classes:
        return empty
    tprs, tnrs, fprs, fnrs = [], [], [], []
    for cls in classes:
        yt = y.eq(cls).to_numpy()
        yp = p.eq(cls).to_numpy()
        tp = float(np.logical_and(yt, yp).sum())
        fn = float(np.logical_and(yt, np.logical_not(yp)).sum())
        fp = float(np.logical_and(np.logical_not(yt), yp).sum())
        tn = float(np.logical_and(np.logical_not(yt), np.logical_not(yp)).sum())
        tprs.append(tp / (tp + fn) if (tp + fn) else 0.0)
        tnrs.append(tn / (tn + fp) if (tn + fp) else 0.0)
        fprs.append(fp / (fp + tn) if (fp + tn) else 0.0)
        fnrs.append(fn / (tp + fn) if (tp + fn) else 0.0)
    tpr = float(np.mean(tprs))
    tnr = float(np.mean(tnrs))
    fpr = float(np.mean(fprs))
    fnr = float(np.mean(fnrs))
    return {
        "tpr": tpr,
        "tnr": tnr,
        "fpr": fpr,
        "fnr": fnr,
        "tn_pct": 100.0 * tnr,
        "fn_pct": 100.0 * fnr,
        "fp_pct": 100.0 * fpr,
        "tp_pct": 100.0 * tpr,
    }


def opal_style_metrics(true: Iterable, pred: Iterable) -> Dict[str, float]:
    """CAMI/OPAL profiling metrics at the current taxid resolution.

    Presence metrics use the sets of classified taxa (unclassified ignored).
    L1 and Bray–Curtis use relative abundances over all classified reads.
    """
    is_special_taxon, normalize_taxon_token = _taxon_helpers()
    y = normalize_series(true)
    p = normalize_series(pred)
    gold = {t for t in pd.unique(y) if not is_special_taxon(t)}
    pred_set = {t for t in pd.unique(p) if not is_special_taxon(t)}
    tp = float(len(gold & pred_set))
    fp = float(len(pred_set - gold))
    fn = float(len(gold - pred_set))
    completeness = tp / (tp + fn) if (tp + fn) else 0.0
    purity = tp / (tp + fp) if (tp + fp) else 0.0
    f1 = (
        0.0
        if completeness + purity == 0
        else 2 * completeness * purity / (completeness + purity)
    )
    union = tp + fp + fn
    jaccard = tp / union if union else 0.0

    y_c = y[~special_mask(y)]
    p_c = p[~special_mask(p)]
    taxa = sorted(gold | pred_set)
    yt = y_c.value_counts(normalize=True) if len(y_c) else pd.Series(dtype=float)
    pt = p_c.value_counts(normalize=True) if len(p_c) else pd.Series(dtype=float)
    l1 = 0.0
    bray_num = 0.0
    bray_den = 0.0
    for t in taxa:
        a = float(yt.get(t, 0.0))
        b = float(pt.get(t, 0.0))
        l1 += abs(a - b)
        bray_num += abs(a - b)
        bray_den += a + b
    bray = bray_num / bray_den if bray_den else 0.0

    out = {
        "opal_tp": tp,
        "opal_fp": fp,
        "opal_fn": fn,
        "completeness": completeness,
        "opal_purity": purity,
        "opal_f1": f1,
        "jaccard": jaccard,
        "l1_norm": l1,
        "bray_curtis": bray,
    }
    out.update(confusion_rates(y, p))
    return out


def write_profiles_for_annotators(
    work: pd.DataFrame,
    annotators: Sequence[str],
    dest_dir: Path,
    true_col: str = "true",
    rank: str = "genus",
) -> Tuple[Optional[Path], List[Tuple[str, Path]]]:
    """Gold + per-annotator CAMI profiles. Returns ``(gold, [(name, path)...])``."""
    if true_col not in work.columns:
        return None, []
    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)
    samples = (
        work["sample"].astype(str)
        if "sample" in work.columns
        else pd.Series(["1"] * len(work), index=work.index)
    )
    gold_parts = []
    for sample, idx in samples.groupby(samples).groups.items():
        gold_parts.append(
            write_cami_profile(
                series_to_counts(work.loc[idx, true_col]),
                dest_dir / f"_gold_{sample}.profile",
                sample_id=str(sample),
                rank=rank,
                total=float(len(idx)),
            ).read_text(encoding="utf-8")
        )
    gold = dest_dir / "goldstandard.profile"
    gold.write_text("\n".join(gold_parts) + "\n", encoding="utf-8")

    profiles: List[Tuple[str, Path]] = []
    seen_names = set()
    for name in annotators:
        if name not in work.columns or name in seen_names:
            continue
        seen_names.add(name)
        chunks = []
        for sample, idx in samples.groupby(samples).groups.items():
            chunks.append(
                write_cami_profile(
                    series_to_counts(work.loc[idx, name]),
                    dest_dir / f"_{name}_{sample}.profile",
                    sample_id=str(sample),
                    rank=rank,
                    total=float(len(idx)),
                ).read_text(encoding="utf-8")
            )
        path = dest_dir / f"{name}.profile"
        path.write_text("\n".join(chunks) + "\n", encoding="utf-8")
        profiles.append((str(name), path))
    return gold, profiles


def complete_opal_metrics(frame: pd.DataFrame) -> pd.DataFrame:
    """Add an F1 row wherever OPAL stored purity and completeness without F1.

    OPAL's ``print_by_tool`` indexes the ``F1 score`` column unconditionally.
    A rank that only has abundance metrics used to abort the whole report
    after ``results.tsv`` was already written.
    """
    needed = {"tool", "rank", "metric", "sample", "value"}
    if not isinstance(frame, pd.DataFrame) or frame.empty or not needed.issubset(frame.columns):
        return frame
    extra = []
    for (tool, rank, sample), part in frame.groupby(["tool", "rank", "sample"], dropna=False):
        metrics = set(part["metric"].astype(str))
        if "F1 score" in metrics or "Purity" not in metrics or "Completeness" not in metrics:
            continue
        try:
            purity = float(part.loc[part["metric"] == "Purity", "value"].iloc[0])
            recall = float(part.loc[part["metric"] == "Completeness", "value"].iloc[0])
        except (TypeError, ValueError):
            continue
        f1 = 0.0 if purity + recall == 0 else 2.0 * purity * recall / (purity + recall)
        extra.append(
            {"tool": tool, "rank": rank, "metric": "F1 score", "sample": sample, "value": f1}
        )
    if not extra:
        return frame
    return pd.concat([frame, pd.DataFrame(extra)], ignore_index=True)


def _dedupe_profiles(
    profiles: Sequence[Tuple[str, Path]],
) -> List[Tuple[str, Path]]:
    """One profile per label. Repeated labels make OPAL's Categorical index fail."""
    seen = set()
    unique: List[Tuple[str, Path]] = []
    for name, path in profiles:
        label = str(name)
        if label in seen:
            continue
        seen.add(label)
        unique.append((label, path))
    return unique


def _opal_score_ready(output_dir: Path) -> bool:
    path = Path(output_dir) / "results.tsv"
    if not path.is_file():
        return False
    try:
        frame = pd.read_csv(path, sep="\t")
    except Exception:
        return False
    return (
        not frame.empty
        and {"metric", "value"}.issubset(frame.columns)
        and frame["metric"].astype(str).ne("").any()
    )


def _install_opal_plot_compat() -> None:
    """Keep OPAL's score tables when a small profile cannot draw a spider chart.

    With three or more tools, a rank that has no values is a length-0 series
    against a length-N axis (``shapes (N,) and (0,)``). That used to exit
    after metrics were saved. Unknown ranks raise ``KeyError`` in the same plot.
    """
    import cami_opal.plots as plots

    if getattr(plots, "_samovar_small_profile", False):
        return
    orig_spider = plots.spider_plot
    orig_plot_all = plots.plot_all

    def spider_plot(metrics, labels, rank_to_metric_to_toolvalues, *args, **kwargs):
        try:
            return orig_spider(metrics, labels, rank_to_metric_to_toolvalues, *args, **kwargs)
        except ValueError as exc:
            if "same first dimension" not in str(exc):
                raise
            logger.info("OPAL spider plot skipped: %s", exc)
            return []

    def plot_all(*args, **kwargs):
        try:
            return orig_plot_all(*args, **kwargs)
        except (KeyError, ValueError) as exc:
            logger.info("OPAL summary plots skipped: %s", exc)
            return []

    plots.spider_plot = spider_plot
    plots.plot_all = plot_all
    plots._samovar_small_profile = True


def _patched_opal_source(script: str) -> str:
    """Tolerate a missing F1 column and repeated category labels in ``opal.py``.

    The installed script is not ours. These are the three small-data crashes:
    ``KeyError: 'F1 score'``, ``Categorical categories must be unique``, and
    the spider-plot shape error handled in ``_install_opal_plot_compat``.
    """
    source = Path(script).read_text(encoding="utf-8")
    source = source.replace(
        "table = table.loc[:, order_columns]",
        "table = table.loc[:, [col for col in order_columns if col in table.columns]]",
    )
    source = source.replace(
        "categories=order_rows, ordered=True",
        "categories=list(dict.fromkeys(order_rows)), ordered=True",
    )
    source = source.replace(
        "categories=ranks_in_metrics, ordered=True",
        "categories=list(dict.fromkeys(ranks_in_metrics)), ordered=True",
    )
    needle = "    logger.info('Saving computed metrics')"
    if needle in source:
        source = source.replace(
            needle,
            "    from samovar.opal import complete_opal_metrics as _complete_opal_metrics\n"
            "    pd_metrics = _complete_opal_metrics(pd_metrics)\n"
            "    logger.info('Saving computed metrics')",
            1,
        )
    return source


def run_opal(
    gold: Path,
    profiles: Sequence[Tuple[str, Path]],
    output_dir: Path,
    rank: str = "genus",
    timeout: int = 180,
) -> Optional[Path]:
    """Invoke OPAL if installed. Returns the output directory on success."""
    profiles = _dedupe_profiles(profiles)
    if not opal_enabled() or not profiles:
        return None
    exe = opal_command()
    if not exe:
        logger.info("OPAL not installed; skipping CAMI HTML (./install.sh OPAL)")
        return None
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    labels = ",".join(name for name, _ in profiles)
    files = [str(path) for _, path in profiles]
    cmd = [
        *exe,
        "-g",
        str(gold),
        "-o",
        str(out),
        "-l",
        labels,
        "-r",
        opal_rank_range(rank),
        "--silent",
        *files,
    ]
    stdout_path, stderr_path = out / "stdout.log", out / "stderr.log"
    try:
        with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
            proc = subprocess.run(cmd, stdout=stdout, stderr=stderr, text=True,
                                  timeout=timeout, check=False)
    except (OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("OPAL did not run: %s (logs: %s, %s)", exc, stdout_path, stderr_path)
        return None
    if proc.returncode != 0:
        if _opal_score_ready(out):
            logger.info(
                "OPAL exited %s after writing results.tsv; score table kept (%s)",
                proc.returncode,
                stderr_path,
            )
            return out
        logger.warning(
            "OPAL exited %s: %s",
            proc.returncode,
            "\n".join(stderr_path.read_text(errors="replace").splitlines()[-30:]) + f" (full log: {stderr_path})",
        )
        return None
    return out


def maybe_run_opal(
    work: pd.DataFrame,
    annotators: Sequence[str],
    output_dir: Path,
    true_col: str = "true",
    rank: str = "genus",
) -> Optional[Path]:
    if true_col not in work.columns:
        return None
    dest = Path(output_dir) / "opal"
    gold, profiles = write_profiles_for_annotators(
        work, annotators, dest / "profiles", true_col=true_col, rank=rank
    )
    if gold is None or not profiles:
        return None
    return run_opal(gold, profiles, dest, rank=rank)


def main():
    """Run the configured OPAL script with process-local matplotlib compatibility."""
    import sys
    import matplotlib
    import matplotlib.cm

    if not hasattr(matplotlib.cm, "get_cmap"):
        matplotlib.cm.get_cmap = matplotlib.colormaps.get_cmap
    script = sys.argv[1]
    sys.argv = sys.argv[1:]
    _install_opal_plot_compat()
    source = _patched_opal_source(script)
    exec(compile(source, script, "exec"), {"__name__": "__main__", "__file__": script})


if __name__ == "__main__":
    main()
