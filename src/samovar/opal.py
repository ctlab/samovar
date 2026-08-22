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
import shutil
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
    from samovar.viz_annotation import is_special_taxon, normalize_taxon_token

    return is_special_taxon, normalize_taxon_token


def cami_rank(rank: Optional[str]) -> str:
    key = "" if rank is None else str(rank).strip().lower()
    mapped = CAMI_RANK_ALIASES.get(key, key or "genus")
    return mapped if mapped in OPAL_RANKS else "genus"


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
    from samovar.paths import load_config, resolve_executable

    cfg = load_config()
    configured = str(cfg.get("opal_path") or "").strip()
    env = os.environ.get("SAMOVAR_OPAL_PATH", "").strip()
    for token in (env, configured, "opal.py", "opal"):
        if not token:
            continue
        resolved = resolve_executable(token, tool_key="opal.py")
        path = (resolved or token).split()[0]
        if path and Path(path).is_file() and os.access(path, os.X_OK):
            return str(Path(path).resolve())
        found = shutil.which(Path(token).name) or shutil.which(token)
        if found:
            return found
    return None


def write_cami_profile(
    counts: Dict[str, float],
    dest: Path,
    sample_id: str = "1",
    rank: str = "genus",
    total: Optional[float] = None,
) -> Path:
    """Write a single-sample CAMI bioboxes taxonomic profile."""
    is_special_taxon, normalize_taxon_token = _taxon_helpers()
    rank = cami_rank(rank)
    classified = float(sum(max(float(v), 0.0) for v in counts.values()))
    denom = float(total) if total is not None and float(total) > 0 else classified
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    sid = "".join(ch if ch.isalnum() or ch in "._" else "_" for ch in str(sample_id)) or "1"
    lines = [
        f"@SampleID:{sid}",
        "@Version:0.9.1",
        f"@Ranks:{rank}",
        "@TaxonomyID:ncbi",
        "@@TAXID\tRANK\tTAXPATH\tPERCENTAGE",
    ]
    for taxid, n in sorted(counts.items(), key=lambda kv: (-float(kv[1]), str(kv[0]))):
        token = normalize_taxon_token(taxid)
        if is_special_taxon(token) or float(n) <= 0:
            continue
        pct = 0.0 if denom <= 0 else 100.0 * float(n) / denom
        lines.append(f"{token}\t{rank}\t{token}\t{pct:.6f}")
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return dest


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
    y = pd.Series(list(true) if not isinstance(true, pd.Series) else true).map(
        normalize_taxon_token
    )
    p = pd.Series(list(pred) if not isinstance(pred, pd.Series) else pred).map(
        normalize_taxon_token
    )
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
    y = pd.Series(list(true) if not isinstance(true, pd.Series) else true).map(
        normalize_taxon_token
    )
    p = pd.Series(list(pred) if not isinstance(pred, pd.Series) else pred).map(
        normalize_taxon_token
    )
    gold = {t for t in y if not is_special_taxon(t)}
    pred_set = {t for t in p if not is_special_taxon(t)}
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

    y_c = y[~y.map(is_special_taxon)]
    p_c = p[~p.map(is_special_taxon)]
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
    for name in annotators:
        if name not in work.columns:
            continue
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


def run_opal(
    gold: Path,
    profiles: Sequence[Tuple[str, Path]],
    output_dir: Path,
    rank: str = "genus",
    timeout: int = 180,
) -> Optional[Path]:
    """Invoke OPAL if installed. Returns the output directory on success."""
    if not opal_enabled() or not profiles:
        return None
    exe = opal_executable()
    if not exe:
        logger.info("OPAL not installed; skipping CAMI HTML (./install.sh OPAL)")
        return None
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    labels = ",".join(name for name, _ in profiles)
    files = [str(path) for _, path in profiles]
    rank_name = cami_rank(rank)
    cmd = [
        exe,
        "-g",
        str(gold),
        "-o",
        str(out),
        "-l",
        labels,
        "-r",
        f"{rank_name},{rank_name}",
        "--silent",
        *files,
    ]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        logger.warning("OPAL did not run: %s", exc)
        return None
    if proc.returncode != 0:
        logger.warning(
            "OPAL exited %s: %s",
            proc.returncode,
            (proc.stderr or proc.stdout or "")[:500],
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
