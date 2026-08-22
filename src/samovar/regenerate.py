"""Regenerate metagenome abundance tables from annotation directories.

Modes mirror the R ``annotation_regenerate.R`` workflow:

- ``preserve`` (default): observed taxID counts, no generative remodelling.
- ``glm``: correlation-aware synthetic communities (Python) or R samovaR when
  ``SAMOVAR_USE_R=1``.
- ``bootstrap``: column bootstrap resampling of observed sample profiles.
- ``vae``: latent-factor generative model (FactorAnalysis on log abundances).
"""

from __future__ import annotations

import os
import re
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd

from samovar.annotation_io import read_annotation_dir

PRESERVE_MODES = frozenset({"preserve", "none", "exact", "raw", "off", "false", ""})
GENERATIVE_MODES = frozenset({"glm", "bootstrap", "vae"})


def normalize_regeneration_mode(mode: Optional[str]) -> str:
    if mode is None:
        return "preserve"
    key = str(mode).strip().lower()
    if key in PRESERVE_MODES:
        return "preserve"
    if key in GENERATIVE_MODES:
        return key
    raise ValueError(
        f"Unknown regeneration_mode={mode!r}. "
        "Use preserve, glm, bootstrap, or vae."
    )


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


def _filter_taxa(matrix: pd.DataFrame, threshold_amount: float, n_reads: Optional[int]) -> pd.DataFrame:
    if matrix.empty:
        return matrix
    mat = matrix.copy()
    if n_reads and n_reads > 0:
        totals = mat.sum(axis=0).replace(0, 1)
        mat = mat / totals * float(n_reads)
    if threshold_amount and threshold_amount > 0:
        if n_reads and n_reads > 0:
            keep = mat.max(axis=1) >= threshold_amount * float(n_reads)
        else:
            keep = mat.max(axis=1) >= threshold_amount
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
        mat = _filter_taxa(mat, threshold_amount, n_reads if rescale else None)
        if rescale and n_reads:
            mat = _scale_columns_to_n_reads(mat, n_reads)
        else:
            mat = mat.round()
        name = _annotator_name(col)
        result[name] = _abundance_table_from_matrix(mat, name)
    return result


def regenerate_bootstrap(
    data: pd.DataFrame,
    n_samples: int,
    n_reads: Optional[int] = None,
    threshold_amount: float = 1e-5,
    seed: Optional[int] = 42,
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
        synth = {}
        for i in range(int(n_samples)):
            picks = rng.integers(0, n_src, size=n_src)
            combined = pd.Series(0.0, index=mat.index)
            for idx in picks:
                combined = combined.add(mat.iloc[:, int(idx)], fill_value=0)
            synth[f"synth_{i + 1}"] = combined
        synth_mat = pd.DataFrame(synth)
        synth_mat = _filter_taxa(synth_mat, threshold_amount, n_reads)
        synth_mat = _scale_columns_to_n_reads(synth_mat, n_reads)
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
        mat = _filter_taxa(mat, threshold_amount, n_reads)
        if mat.shape[0] < 2 or mat.shape[1] < 2:
            name = _annotator_name(col)
            result[name] = _abundance_table_from_matrix(mat.round(), name)
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
        for i in range(int(n_samples)):
            z = rng.normal(0, 1, size=n_comp) * latent_std
            decoded = fa.mean_ + np.dot(z.reshape(1, -1), fa.components_).flatten()
            decoded = np.expm1(np.maximum(decoded, 0))
            synth_cols[f"synth_{i + 1}"] = decoded
        synth_mat = pd.DataFrame(synth_cols, index=mat.index)
        synth_mat = _scale_columns_to_n_reads(synth_mat, n_reads)
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
        mat = _filter_taxa(mat, threshold_amount, n_reads)
        if mat.empty:
            continue
        log_mat = np.log1p(mat.to_numpy(dtype=float))
        corr = np.corrcoef(log_mat)
        np.fill_diagonal(corr, 1.0)
        synth = {}
        src_cols = list(range(mat.shape[1]))
        for i in range(int(n_samples)):
            base_idx = int(rng.choice(src_cols))
            profile = log_mat[:, base_idx].copy()
            for j in range(profile.shape[0]):
                partners = corr[j] * rng.normal(0, noise_scale, size=profile.shape[0])
                profile[j] += float(np.sum(partners)) / max(profile.shape[0], 1)
            decoded = np.expm1(np.maximum(profile, 0))
            synth[f"synth_{i + 1}"] = decoded
        synth_mat = pd.DataFrame(synth, index=mat.index)
        synth_mat = _scale_columns_to_n_reads(synth_mat, n_reads)
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
    mode = normalize_regeneration_mode(cfg.get("regeneration_mode", "preserve"))
    n_samples = int(cfg.get("N", 10))
    n_reads = cfg.get("N_reads")
    if n_reads is not None:
        n_reads = int(n_reads)
    threshold = float(cfg.get("threshold_amount", cfg.get("treshhold_amount", 1e-5)))
    seed = cfg.get("seed", 42)
    rescale = bool(cfg.get("rescale_abundance", False))
    latent_dim = int(cfg.get("vae_latent_dim", 4))

    if data is None:
        data = read_annotation_dir(annotation_dir)

    if mode == "preserve":
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
        )
    elif mode == "vae":
        tables = regenerate_vae(
            data,
            n_samples=n_samples,
            n_reads=n_reads,
            threshold_amount=threshold,
            latent_dim=latent_dim,
            seed=seed,
        )
    else:  # glm
        tables = regenerate_glm_python(
            data,
            n_samples=n_samples,
            n_reads=n_reads,
            threshold_amount=threshold,
            seed=seed,
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
        "N": 10,
        "N_reads": 1000,
        "regeneration_mode": "preserve",
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
