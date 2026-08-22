"""Regenerate metagenome abundance tables from annotation directories.

Modes:

- ``direct`` (default; aliases: preserve, exact, raw): observed taxID counts,
  same sample names, no generative remodelling.
- ``bootstrap``: column bootstrap of observed sample profiles.
- ``vae``: latent-factor generative model (FactorAnalysis on log abundances).
- ``glm``: correlation-aware synthetic communities (Python).
- ``samovar``: optional R regenerator (not part of the Python install). Looked
  up via ``SAMOVAR_R_REGENERATE`` / config ``annotation_regenerate_r``.
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

from samovar.annotation_io import read_annotation_dir

DIRECT_MODES = frozenset(
    {"direct", "preserve", "none", "exact", "raw", "off", "false", ""}
)
GENERATIVE_MODES = frozenset({"glm", "bootstrap", "vae"})
SAMOVAR_R_MODES = frozenset({"samovar", "r", "boil"})


def normalize_regeneration_mode(mode: Optional[str]) -> str:
    if mode is None:
        return "direct"
    key = str(mode).strip().lower()
    if key in DIRECT_MODES:
        return "direct"
    if key in GENERATIVE_MODES:
        return key
    if key in SAMOVAR_R_MODES:
        return "samovar"
    raise ValueError(
        f"Unknown regeneration_mode={mode!r}. "
        "Use direct, bootstrap, vae, glm, or samovar."
    )


def is_direct_mode(mode: Optional[str]) -> bool:
    return normalize_regeneration_mode(mode) == "direct"


def _taxid_columns(df: pd.DataFrame) -> List[str]:
    cols = []
    for col in df.columns:
        name = str(col)
        if name in {"seq", "length", "sample", "true"}:
            continue
        if "confidence" in name.lower() or name.endswith("_conf"):
            continue
        if name.startswith("taxID") or name == "true":
            cols.append(col)
    return cols


def _annotator_name(col: str) -> str:
    name = re.sub(r"^taxID_", "", str(col))
    name = re.sub(r"_[0-9]+$", "", name)
    return name


def _count_matrix(
    data: pd.DataFrame,
    tax_col: str,
    sample_col: str = "sample",
) -> pd.DataFrame:
    """Species x sample count matrix (taxid rows, sample columns)."""
    work = data[[tax_col, sample_col]].copy()
    work[tax_col] = work[tax_col].astype(str)
    counts = work.groupby([tax_col, sample_col], dropna=False).size().unstack(fill_value=0)
    counts.index = counts.index.astype(str)
    return counts


def _filter_taxa(
    matrix: pd.DataFrame,
    threshold_amount: float,
    n_reads: Optional[int] = None,
    rescale: bool = False,
) -> pd.DataFrame:
    if matrix.empty:
        return matrix
    mat = matrix.copy()
    if rescale and n_reads and n_reads > 0:
        totals = mat.sum(axis=0).replace(0, 1)
        mat = mat / totals * float(n_reads)
        if threshold_amount and threshold_amount > 0:
            keep = mat.max(axis=1) >= threshold_amount * float(n_reads)
            mat = mat.loc[keep]
        return mat
    if threshold_amount and threshold_amount > 0:
        totals = mat.sum(axis=0).replace(0, 1)
        fracs = mat / totals
        keep = fracs.max(axis=1) >= float(threshold_amount)
        if not bool(keep.any()):
            keep = mat.max(axis=1) > 0
        mat = mat.loc[keep]
    return mat


def _scale_columns_to_n_reads(matrix: pd.DataFrame, n_reads: Optional[int]) -> pd.DataFrame:
    if not n_reads or matrix.empty:
        return matrix.round()
    out = matrix.copy()
    target = int(n_reads)
    for col in out.columns:
        total = out[col].sum()
        if total <= 0:
            continue
        scaled = (out[col] / total * float(target)).round()
        diff = target - int(scaled.sum())
        if diff != 0 and scaled.sum() > 0:
            idx = scaled.idxmax()
            scaled.loc[idx] = scaled.loc[idx] + diff
        out[col] = scaled
    return out.round()


def _apply_abundance_scale(
    matrix: pd.DataFrame,
    n_reads: Optional[int],
    rescale: bool,
) -> pd.DataFrame:
    if rescale:
        return _scale_columns_to_n_reads(matrix, n_reads)
    return matrix.round()


def _abundance_table_from_matrix(matrix: pd.DataFrame, annotator_name: str) -> pd.DataFrame:
    out = matrix.copy()
    out = out.reset_index()
    out.columns = ["taxid"] + [f"N_{c}" for c in out.columns[1:]]
    # Drop all-zero taxa
    n_cols = [c for c in out.columns if c.startswith("N_")]
    if n_cols:
        out = out[out[n_cols].sum(axis=1) > 0]
    return out


def regenerate_preserve(
    data: pd.DataFrame,
    n_reads: Optional[int] = None,
    rescale: bool = False,
    threshold_amount: float = 0.0,
) -> Dict[str, pd.DataFrame]:
    """Count taxIDs per sample; default keeps observed abundances."""
    result: Dict[str, pd.DataFrame] = {}
    sample_col = "sample" if "sample" in data.columns else None
    if sample_col is None:
        data = data.copy()
        data["sample"] = "1"
        sample_col = "sample"
    for col in _taxid_columns(data):
        mat = _count_matrix(data, col, sample_col)
        mat = _filter_taxa(mat, threshold_amount, n_reads, rescale=rescale)
        mat = _apply_abundance_scale(mat, n_reads, rescale)
        name = _annotator_name(col)
        result[name] = _abundance_table_from_matrix(mat, name)
    return result


def regenerate_bootstrap(
    data: pd.DataFrame,
    n_samples: int,
    n_reads: Optional[int] = None,
    threshold_amount: float = 1e-5,
    seed: Optional[int] = 42,
    rescale: bool = False,
) -> Dict[str, pd.DataFrame]:
    """Bootstrap resample observed sample columns to synthesize new profiles."""
    rng = np.random.default_rng(seed)
    result: Dict[str, pd.DataFrame] = {}
    sample_col = "sample" if "sample" in data.columns else "sample"
    work = data.copy()
    if sample_col not in work.columns:
        work[sample_col] = "1"
    for col in _taxid_columns(work):
        mat = _count_matrix(work, col, sample_col)
        if mat.shape[1] == 0:
            continue
        n_src = mat.shape[1]
        names = synthetic_sample_names(list(mat.columns), n_samples)
        synth = {}
        for name in names:
            picks = rng.integers(0, n_src, size=n_src)
            combined = pd.Series(0.0, index=mat.index)
            for idx in picks:
                combined = combined.add(mat.iloc[:, int(idx)], fill_value=0)
            synth[name] = combined
        synth_mat = pd.DataFrame(synth)
        synth_mat = _filter_taxa(synth_mat, threshold_amount, n_reads, rescale=rescale)
        synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
        name = _annotator_name(col)
        result[name] = _abundance_table_from_matrix(synth_mat, name)
    return result


def regenerate_vae(
    data: pd.DataFrame,
    n_samples: int,
    n_reads: Optional[int] = 1000,
    threshold_amount: float = 1e-5,
    latent_dim: int = 4,
    seed: Optional[int] = 42,
    rescale: bool = False,
) -> Dict[str, pd.DataFrame]:
    """Sample synthetic profiles with a latent linear generative model (FA)."""
    from sklearn.decomposition import FactorAnalysis

    rng = np.random.default_rng(seed)
    result: Dict[str, pd.DataFrame] = {}
    sample_col = "sample" if "sample" in data.columns else "sample"
    work = data.copy()
    if sample_col not in work.columns:
        work[sample_col] = "1"
    for col in _taxid_columns(work):
        mat = _count_matrix(work, col, sample_col)
        mat = _filter_taxa(mat, threshold_amount, n_reads, rescale=rescale)
        names = synthetic_sample_names(list(mat.columns), n_samples)
        if mat.shape[0] < 2 or mat.shape[1] < 2:
            synth = {}
            cols = list(mat.columns) or ["1"]
            for i, name in enumerate(names):
                src = cols[i % len(cols)]
                synth[name] = mat[src] if src in mat.columns else mat.iloc[:, 0]
            synth_mat = pd.DataFrame(synth, index=mat.index)
            synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
            name = _annotator_name(col)
            result[name] = _abundance_table_from_matrix(synth_mat, name)
            continue
        log_mat = np.log1p(mat.to_numpy(dtype=float))
        X = log_mat.T  # sklearn: samples x features (taxa)
        n_comp = min(int(latent_dim), X.shape[0] - 1, X.shape[1] - 1)
        n_comp = max(1, n_comp)
        fa = FactorAnalysis(n_components=n_comp, random_state=int(seed or 0))
        fa.fit(X)
        latent = fa.transform(X)
        latent_std = np.std(latent, axis=0)
        latent_std[latent_std == 0] = 1.0
        synth_cols = {}
        for name in names:
            z = rng.normal(0, 1, size=n_comp) * latent_std
            decoded = fa.mean_ + np.dot(z.reshape(1, -1), fa.components_).flatten()
            decoded = np.expm1(np.maximum(decoded, 0))
            synth_cols[name] = decoded
        synth_mat = pd.DataFrame(synth_cols, index=mat.index)
        synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
        name = _annotator_name(col)
        result[name] = _abundance_table_from_matrix(synth_mat, name)
    return result


def regenerate_glm_python(
    data: pd.DataFrame,
    n_samples: int,
    n_reads: Optional[int] = 1000,
    threshold_amount: float = 1e-5,
    seed: Optional[int] = 42,
    noise_scale: float = 0.15,
    rescale: bool = False,
) -> Dict[str, pd.DataFrame]:
    """Correlation-perturbed resampling (Python analog of glm-based generation)."""
    rng = np.random.default_rng(seed)
    result: Dict[str, pd.DataFrame] = {}
    sample_col = "sample" if "sample" in data.columns else "sample"
    work = data.copy()
    if sample_col not in work.columns:
        work[sample_col] = "1"
    for col in _taxid_columns(work):
        mat = _count_matrix(work, col, sample_col)
        mat = _filter_taxa(mat, threshold_amount, n_reads, rescale=rescale)
        if mat.empty:
            continue
        names = synthetic_sample_names(list(mat.columns), n_samples)
        log_mat = np.log1p(mat.to_numpy(dtype=float))
        corr = _correlation_matrix(log_mat)
        synth = {}
        src_cols = list(range(mat.shape[1]))
        if not src_cols:
            continue
        n_taxa = int(log_mat.shape[0])
        for name in names:
            base_idx = int(rng.choice(src_cols))
            profile = log_mat[:, base_idx].copy()
            if n_taxa >= 1 and corr.shape[0] == n_taxa:
                noise = rng.normal(0, noise_scale, size=n_taxa)
                for j in range(n_taxa):
                    partners = corr[j] * noise
                    profile[j] += float(np.sum(partners)) / max(n_taxa, 1)
            decoded = np.expm1(np.maximum(profile, 0))
            synth[name] = decoded
        synth_mat = pd.DataFrame(synth, index=mat.index)
        synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
        name = _annotator_name(col)
        result[name] = _abundance_table_from_matrix(synth_mat, name)
    return result


def regenerate_annotation_tables(
    annotation_dir: Union[str, Path],
    output_dir: Union[str, Path],
    config: Optional[Dict[str, Any]] = None,
    data: Optional[pd.DataFrame] = None,
) -> Dict[str, pd.DataFrame]:
    """Regenerate per-annotator abundance CSVs from an annotation directory."""
    cfg = dict(config or {})
    mode = normalize_regeneration_mode(cfg.get("regeneration_mode", "direct"))
    n_reads = cfg.get("N_reads")
    if n_reads is not None:
        n_reads = int(n_reads)
    threshold = float(cfg.get("threshold_amount", cfg.get("treshhold_amount", 1e-5)))
    seed = cfg.get("seed", 42)
    rescale = bool(cfg.get("rescale_abundance", False))
    latent_dim = int(cfg.get("vae_latent_dim", 4))

    if data is None:
        data = read_annotation_dir(annotation_dir)
    n_samples = _n_samples_or_observed(data, cfg.get("N"))

    if mode == "direct":
        tables = regenerate_preserve(
            data,
            n_reads=n_reads,
            rescale=rescale,
            threshold_amount=threshold,
        )
    elif mode == "bootstrap":
        tables = regenerate_bootstrap(
            data,
            n_samples=n_samples,
            n_reads=n_reads,
            threshold_amount=threshold,
            seed=seed,
            rescale=rescale,
        )
    elif mode == "vae":
        tables = regenerate_vae(
            data,
            n_samples=n_samples,
            n_reads=n_reads,
            threshold_amount=threshold,
            latent_dim=latent_dim,
            seed=seed,
            rescale=rescale,
        )
    elif mode == "glm":
        tables = regenerate_glm_python(
            data,
            n_samples=n_samples,
            n_reads=n_reads,
            threshold_amount=threshold,
            seed=seed,
            rescale=rescale,
        )
    else:
        raise ValueError(
            "regeneration_mode='samovar' uses the optional R regenerator. "
            "Call samovar_annotation_regenerate() after installing it, or set "
            "SAMOVAR_R_REGENERATE to annotation_regenerate.R."
        )

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name))
        table.to_csv(out / f"{safe}.csv", index=False)
    return tables


def write_samovar_config_defaults(config: Dict[str, Any]) -> Dict[str, Any]:
    """Merge SamovaR-style defaults used by R ``annotation_regenerate.R``."""
    defaults = {
        "threshold_amount": 1e-5,
        "treshhold_amount": 1e-5,
        "plot_log": False,
        "min_cluster_size": 2,
        "max_cluster_size": 100,
        "N_reads": 1000,
        "regeneration_mode": "direct",
        "seed": 42,
        "vae_latent_dim": 4,
        "rescale_abundance": False,
    }
    merged = {**defaults, **config}
    # keep R spelling alias
    if "treshhold_amount" in merged and "threshold_amount" not in config:
        merged["threshold_amount"] = merged["treshhold_amount"]
    return merged


def sample_names_from_abundance_columns(
    n_columns: List[str],
    sample_names_hint: Optional[List[str]] = None,
) -> List[str]:
    """Map ``N_*`` abundance columns to pipeline sample names."""
    if sample_names_hint and len(sample_names_hint) == len(n_columns):
        return list(sample_names_hint)
    names: List[str] = []
    for col in n_columns:
        col_s = str(col)
        if col_s.startswith("N_"):
            names.append(col_s[2:])
        else:
            names.append(col_s)
    return names


def synthetic_sample_names(
    original: Sequence[str],
    n_samples: Optional[int] = None,
) -> List[str]:
    """Stable sample names for generative regeneration.

    Names are derived from the observed samples. ``N`` only changes how many
    profiles are emitted (truncate, or extra rounds as ``{name}_r2``), never
    switches the scheme to ``synth_1..N``.
    """
    orig: List[str] = []
    seen = set()
    for name in original:
        text = str(name).strip()
        if not text or text in seen:
            continue
        orig.append(text)
        seen.add(text)
    if not orig:
        orig = ["1"]
    if n_samples is None or int(n_samples) <= 0:
        n = len(orig)
    else:
        n = int(n_samples)
    if n <= len(orig):
        return orig[:n]
    names = list(orig)
    round_idx = 2
    i = 0
    while len(names) < n:
        names.append(f"{orig[i % len(orig)]}_r{round_idx}")
        i += 1
        if i % len(orig) == 0:
            round_idx += 1
    return names


def _n_samples_or_observed(data: pd.DataFrame, n_samples: Optional[int]) -> int:
    observed = 1
    if "sample" in data.columns and len(data.index):
        observed = max(int(pd.Series(data["sample"].astype(str)).nunique()), 1)
    if n_samples is None or int(n_samples) <= 0:
        return observed
    return int(n_samples)


def _correlation_matrix(log_mat: np.ndarray) -> np.ndarray:
    """Row-wise correlation; always finite with 1s on the diagonal."""
    n_taxa = int(log_mat.shape[0]) if log_mat.ndim == 2 else 0
    if n_taxa <= 0:
        return np.zeros((0, 0), dtype=float)
    if n_taxa == 1 or log_mat.shape[1] < 2:
        return np.ones((n_taxa, n_taxa), dtype=float)
    with np.errstate(invalid="ignore", divide="ignore"):
        corr = np.corrcoef(log_mat)
    corr = np.atleast_2d(np.asarray(corr, dtype=float))
    if corr.shape != (n_taxa, n_taxa):
        return np.eye(n_taxa, dtype=float)
    corr = np.nan_to_num(corr, nan=0.0, posinf=0.0, neginf=0.0)
    np.fill_diagonal(corr, 1.0)
    return corr
