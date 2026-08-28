"""Regenerate metagenome abundance tables from annotation directories.

Modes:

- ``direct`` (default; aliases: preserve, exact, raw): observed taxID counts,
  same sample names, no generative remodelling.
- ``bootstrap``: column bootstrap of observed sample profiles.
- ``vae``: latent-factor generative model (FactorAnalysis on log abundances).
- ``glm``: correlation-aware synthetic communities (Python).
- ``camisim-table`` (aliases camisim, cami): CAMISIM log-normal community
  design on the observed taxIDs (same engine as ``samovar generate --simulator camisim --camisim-mode table``).
- ``samovar``: optional R regenerator (not part of the Python install). Looked
  up via ``SAMOVAR_R_REGENERATE`` / config ``annotation_regenerate_r``.
- any other name: custom ``table_reads_generator`` registered with
  ``samovar tools import --type table`` (see ``table_regenerators``).
"""

from __future__ import annotations

import os
import re
import zlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from samovar.annotation_io import read_annotation_dir

DIRECT_MODES = frozenset(
    {"direct", "preserve", "none", "exact", "raw", "off", "false", ""}
)
GENERATIVE_MODES = frozenset({"glm", "bootstrap", "vae"})
SAMOVAR_R_MODES = frozenset({"samovar", "r", "boil"})
CAMISIM_TABLE_MODES = frozenset(
    {"camisim", "camisim-table", "camisim_table", "cami"}
)


def resolve_regeneration_mode(mode: Optional[str]) -> Tuple[str, str]:
    """Return ``("builtin", canonical)`` or ``("custom", name)``.

    Built-ins: ``direct`` (aliases preserve/exact/…), ``bootstrap``, ``vae``,
    ``glm``, ``samovar`` (aliases r/boil). Any other non-empty string is a
    custom ``table_reads_generator`` name (validated against the install
    config at prepare / run).
    """
    if mode is None:
        return "builtin", "direct"
    key = str(mode).strip()
    if not key:
        return "builtin", "direct"
    low = key.lower()
    if low in DIRECT_MODES:
        return "builtin", "direct"
    if low in GENERATIVE_MODES:
        return "builtin", low
    if low in SAMOVAR_R_MODES:
        return "builtin", "samovar"
    if low in CAMISIM_TABLE_MODES:
        return "builtin", "camisim-table"
    return "custom", key


def normalize_regeneration_mode(mode: Optional[str]) -> str:
    """Built-in modes only. Custom names must go through ``resolve_regeneration_mode``."""
    kind, name = resolve_regeneration_mode(mode)
    if kind == "builtin":
        return name
    raise ValueError(
        f"Unknown regeneration_mode={mode!r}. "
        "Use direct, bootstrap, vae, glm, camisim-table, samovar, or an imported "
        "table_reads_generator name (`samovar tools import --type table`)."
    )


def is_direct_mode(mode: Optional[str]) -> bool:
    kind, name = resolve_regeneration_mode(mode)
    return kind == "builtin" and name == "direct"


def coerce_seed(seed: Any, default: int = 42) -> int:
    """Stable non-negative seed for numpy / sklearn / ISS."""
    if seed is None or seed is False:
        return int(default)
    try:
        return int(seed)
    except (TypeError, ValueError):
        return int(default)


def rng_for(seed: Any, *parts: Any) -> np.random.Generator:
    """Independent reproducible stream per (seed, annotator, sample, ...)."""
    ints = [coerce_seed(seed)]
    for part in parts:
        if isinstance(part, (int, np.integer)):
            ints.append(int(part) & 0xFFFFFFFF)
        else:
            ints.append(zlib.crc32(str(part).encode("utf-8")) & 0xFFFFFFFF)
    return np.random.default_rng(np.random.SeedSequence(ints))


def _noisy_copy_column(
    counts: np.ndarray,
    rng: np.random.Generator,
    error_scale: float = 0.15,
) -> np.ndarray:
    """Dirichlet-multinomial resample of one abundance column (not a clone)."""
    counts = np.maximum(np.asarray(counts, dtype=float), 0.0)
    total = int(round(float(counts.sum())))
    if total <= 0:
        return np.zeros_like(counts, dtype=float)
    concentration = max(1.0, 1.0 / max(float(error_scale), 1e-6))
    alpha = counts * (concentration / float(total)) + 0.05
    p = rng.dirichlet(alpha)
    p = p / p.sum()
    return rng.multinomial(total, p).astype(float)


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
    error_scale: float = 0.15,
) -> Dict[str, pd.DataFrame]:
    """Noisy copies of observed samples (Dirichlet-multinomial resampling).

    Each synthetic column is a resample of a real sample with sampling error,
    not a mixture of all columns. Named samples prefer the matching source.
    A fixed seed is reproducible but does not clone every column: each sample
    name gets its own RNG stream.
    """
    seed = coerce_seed(seed)
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
        src_names = [str(c) for c in mat.columns]
        src_index = {name: i for i, name in enumerate(src_names)}
        names = synthetic_sample_names(src_names, n_samples)
        values = mat.to_numpy(dtype=float)
        synth = {}
        annotator = _annotator_name(col)
        pick_rng = rng_for(seed, annotator, "source-pick")
        for name in names:
            rng = rng_for(seed, annotator, name)
            if name in src_index:
                src = src_index[name]
            else:
                src = int(pick_rng.integers(0, n_src))
            synth[name] = pd.Series(
                _noisy_copy_column(values[:, src], rng, error_scale=error_scale),
                index=mat.index,
            )
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

    seed = coerce_seed(seed)
    result: Dict[str, pd.DataFrame] = {}
    sample_col = "sample" if "sample" in data.columns else "sample"
    work = data.copy()
    if sample_col not in work.columns:
        work[sample_col] = "1"
    for col in _taxid_columns(work):
        mat = _count_matrix(work, col, sample_col)
        mat = _filter_taxa(mat, threshold_amount, n_reads, rescale=rescale)
        names = synthetic_sample_names(list(mat.columns), n_samples)
        annotator = _annotator_name(col)
        if mat.shape[0] < 2 or mat.shape[1] < 2:
            synth = {}
            cols = list(mat.columns) or ["1"]
            for i, name in enumerate(names):
                src = cols[i % len(cols)]
                src_vec = mat[src] if src in mat.columns else mat.iloc[:, 0]
                synth[name] = _noisy_copy_column(
                    src_vec.to_numpy(dtype=float),
                    rng_for(seed, annotator, name),
                )
            synth_mat = pd.DataFrame(synth, index=mat.index)
            synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
            result[annotator] = _abundance_table_from_matrix(synth_mat, annotator)
            continue
        log_mat = np.log1p(mat.to_numpy(dtype=float))
        X = log_mat.T  # sklearn: samples x features (taxa)
        n_comp = min(int(latent_dim), X.shape[0] - 1, X.shape[1] - 1)
        n_comp = max(1, n_comp)
        fa = FactorAnalysis(n_components=n_comp, random_state=seed)
        fa.fit(X)
        latent = fa.transform(X)
        latent_std = np.std(latent, axis=0)
        latent_std[latent_std == 0] = 1.0
        synth_cols = {}
        for name in names:
            rng = rng_for(seed, annotator, name)
            z = rng.normal(0, 1, size=n_comp) * latent_std
            decoded = fa.mean_ + np.dot(z.reshape(1, -1), fa.components_).flatten()
            decoded = np.expm1(np.maximum(decoded, 0))
            synth_cols[name] = decoded
        synth_mat = pd.DataFrame(synth_cols, index=mat.index)
        synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
        result[annotator] = _abundance_table_from_matrix(synth_mat, annotator)
    return result


def regenerate_glm_python(
    data: pd.DataFrame,
    n_samples: int,
    n_reads: Optional[int] = 1000,
    threshold_amount: float = 1e-5,
    seed: Optional[int] = 42,
    noise_scale: float = 0.15,
    rescale: bool = False,
    min_cluster_size: int = 2,
    max_cluster_size: int = 100,
) -> Dict[str, pd.DataFrame]:
    """Cluster + Gaussian GLM walk matching ``samovar_boil`` / ``concotion_pour``.

    Taxa are clustered, ``y ~ x`` Gaussian GLMs are fit on co-occurring samples
    (oriented co-occurrence, Rsq edge weights), then each synthetic community
    is generated by walking that graph from a random init taxon.
    """
    seed = coerce_seed(seed)
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
        values = np.maximum(mat.to_numpy(dtype=float), 0.0)
        synth_cols = {}
        annotator = _annotator_name(col)
        for name in names:
            profile = _glm_boil_one(
                values,
                rng_for(seed, annotator, name),
                min_cluster_size=min_cluster_size,
                max_cluster_size=max_cluster_size,
                noise_scale=noise_scale,
            )
            synth_cols[name] = profile
        synth_mat = pd.DataFrame(synth_cols, index=mat.index)
        synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
        result[annotator] = _abundance_table_from_matrix(synth_mat, annotator)
    return result


def regenerate_camisim(
    data: pd.DataFrame,
    n_samples: int,
    n_reads: Optional[int] = 1000,
    threshold_amount: float = 1e-5,
    seed: Optional[int] = 42,
    rescale: bool = False,
    distribution: str = "differential",
    log_mu: float = 1.0,
    log_sigma: float = 2.0,
    gauss_mu: float = 1.0,
    gauss_sigma: float = 1.0,
) -> Dict[str, pd.DataFrame]:
    """CAMISIM-style log-normal abundances for taxIDs seen in the annotation table."""
    from samovar.camisim import design_communities

    seed = coerce_seed(seed)
    depth = int(n_reads) if n_reads not in (None, "", 0) else 1000
    result: Dict[str, pd.DataFrame] = {}
    work = data.copy()
    sample_col = "sample"
    if sample_col not in work.columns:
        work[sample_col] = "1"
    for col in _taxid_columns(work):
        mat = _count_matrix(work, col, sample_col)
        mat = _filter_taxa(mat, threshold_amount, n_reads, rescale=False)
        if mat.empty:
            continue
        taxids = [str(t) for t in mat.index if str(t) not in {"0", "nan", "None", ""}]
        if not taxids:
            continue
        names = synthetic_sample_names(list(mat.columns), n_samples)
        annotator = _annotator_name(col)
        communities = design_communities(
            taxids,
            n_samples=len(names),
            seed=seed + (zlib.crc32(annotator.encode()) & 0xFFFF),
            mode=str(distribution or "differential"),
            log_mu=float(log_mu),
            log_sigma=float(log_sigma),
            gauss_mu=float(gauss_mu),
            gauss_sigma=float(gauss_sigma),
        )
        synth_cols = {}
        for name, frac in zip(names, communities):
            vec = np.array(
                [max(float(frac.get(tid, 0.0)), 0.0) for tid in taxids], dtype=float
            )
            total = float(vec.sum())
            if total <= 0:
                vec = np.ones(len(taxids), dtype=float)
                total = float(vec.sum())
            counts = np.round(vec / total * float(depth))
            if counts.sum() <= 0 and len(counts):
                counts[int(np.argmax(vec))] = float(depth)
            synth_cols[name] = counts
        synth_mat = pd.DataFrame(synth_cols, index=taxids)
        synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
        result[annotator] = _abundance_table_from_matrix(synth_mat, annotator)
    return result


def _write_abundance_tables(
    output_dir: Union[str, Path],
    tables: Dict[str, pd.DataFrame],
) -> None:
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    for name, table in tables.items():
        safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name))
        table.to_csv(out / f"{safe}.csv", index=False)


def _regenerate_one_mode(
    annotation_dir: Union[str, Path],
    output_dir: Union[str, Path],
    cfg: Dict[str, Any],
    data: pd.DataFrame,
) -> Dict[str, pd.DataFrame]:
    from samovar.table_regenerators import (
        as_annotation,
        attach_regenerator_flags,
        get_table_regenerator,
        load_samples_metadata,
    )

    mode = cfg.get("table_reads_generator") or cfg.get("regeneration_mode", "direct")
    local = dict(cfg)
    local["annotation_dir"] = str(annotation_dir)
    local["output_dir"] = str(output_dir)
    local = attach_regenerator_flags(mode, local)
    regenerator = get_table_regenerator(mode)
    tables = regenerator.run(
        as_annotation(data, annotation_dir),
        load_samples_metadata(local),
        local,
    )
    _write_abundance_tables(output_dir, tables)
    return tables


def regenerate_annotation_tables(
    annotation_dir: Union[str, Path],
    output_dir: Union[str, Path],
    config: Optional[Dict[str, Any]] = None,
    data: Optional[pd.DataFrame] = None,
) -> Dict[str, pd.DataFrame]:
    """Regenerate per-annotator abundance CSVs from an annotation directory.

    If ``table_reads_generators`` (or a list ``table_reads_generator``) names
    more than one method, each candidate is scored with ``table_score``
    (default ``shannon_ks``) and only the winner is written to ``output_dir``.
    """
    from samovar.table_regenerators import canonical_regeneration_modes
    from samovar.table_scorers import (
        pick_best_table_method,
        score_generated_tables,
        write_table_selection,
    )

    cfg = dict(config or {})
    if data is None:
        data = read_annotation_dir(annotation_dir)
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    raw_modes = (
        cfg.get("table_reads_generators")
        if cfg.get("table_reads_generators") not in (None, "", False, [])
        else cfg.get("table_reads_generator") or cfg.get("regeneration_mode", "direct")
    )
    modes = canonical_regeneration_modes(raw_modes) or ["direct"]

    if len(modes) == 1:
        one = dict(cfg)
        one["table_reads_generator"] = modes[0]
        one["regeneration_mode"] = modes[0]
        return _regenerate_one_mode(annotation_dir, out, one, data)

    scorer = cfg.get("table_score") or cfg.get("table_reads_scorer") or "shannon_ks"
    candidates_root = out / ".table_candidates"
    rows: List[Dict[str, Any]] = []
    tables_by_mode: Dict[str, Dict[str, pd.DataFrame]] = {}
    for mode in modes:
        cand = candidates_root / re.sub(r"[^A-Za-z0-9._-]+", "_", mode)
        one = dict(cfg)
        one["table_reads_generator"] = mode
        one["regeneration_mode"] = mode
        one["table_reads_generators"] = [mode]
        try:
            tables = _regenerate_one_mode(annotation_dir, cand, one, data)
            score = score_generated_tables(data, tables, scorer)
            score["ok"] = True
        except Exception as exc:
            tables = {}
            score = {
                "scorer": str(scorer),
                "ks_statistic": 1.0,
                "pvalue": 0.0,
                "rank_value": 1.0,
                "ok": False,
                "error": str(exc),
            }
        score["mode"] = mode
        rows.append(score)
        tables_by_mode[mode] = tables

    winner = pick_best_table_method(rows) or modes[0]
    winning_tables = tables_by_mode.get(winner) or {}
    if not winning_tables:
        raise RuntimeError(
            "No table_reads_generator produced abundance tables. "
            f"Tried: {', '.join(modes)}"
        )
    _write_abundance_tables(out, winning_tables)
    payload = {
        "scorer": rows[0].get("scorer", scorer) if rows else scorer,
        "winner": winner,
        "candidates": rows,
    }
    write_table_selection(out / "table_selection.json", payload)
    write_table_selection(candidates_root / "table_selection.json", payload)
    return winning_tables


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
        "bootstrap_error_scale": 0.15,
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


def _occurrence(vec: np.ndarray, min_value: float = 0.0) -> np.ndarray:
    return np.asarray(vec, dtype=float) > float(min_value)


def _oriented_prob(x: np.ndarray, y: np.ndarray, min_value: float = 0.0) -> float:
    """P(Y|X) as in R ``prob_calc_general(..., probability_calculation='oriented')``."""
    x_on = _occurrence(x, min_value)
    n_x = int(x_on.sum())
    if n_x == 0:
        return 0.0
    return float(_occurrence(y[x_on], min_value).sum() / n_x)


def _fit_gaussian_glm(x: np.ndarray, y: np.ndarray, min_value: float = 0.0):
    """Gaussian GLM ``y ~ x`` on co-occurring samples (R ``glm_calc_general``)."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    mask = _occurrence(x, min_value) & _occurrence(y, min_value)
    n = int(mask.sum())
    if n == 0:
        return ("zero",)
    if n == 1:
        xv = float(x[mask][0])
        yv = float(y[mask][0])
        return ("ratio", xv / yv if yv != 0 else 0.0)
    A = np.column_stack([np.ones(n), x[mask]])
    coef, *_ = np.linalg.lstsq(A, y[mask], rcond=None)
    yhat = A @ coef
    ss_res = float(np.sum((y[mask] - yhat) ** 2))
    ymean = float(np.mean(y[mask]))
    ss_tot = float(np.sum((y[mask] - ymean) ** 2))
    rsq = 1e-5 if ss_tot <= 0 else max(1e-5, 1.0 - ss_res / ss_tot)
    return ("ols", float(coef[0]), float(coef[1]), float(rsq))


def _predict_gaussian_glm(fit, xval: float) -> float:
    """R ``get_newval_general`` for mode glm."""
    if fit[0] == "zero":
        return 0.0
    if fit[0] == "ratio":
        return float(fit[1]) * float(xval)
    return max(0.0, float(fit[1]) + float(fit[2]) * float(xval))


def _glm_rsq(fit) -> float:
    if fit[0] == "ols":
        return float(fit[3])
    return 1e-5


def _cluster_taxa(values: np.ndarray, min_cluster_size: int, max_cluster_size: int) -> np.ndarray:
    n_taxa = int(values.shape[0])
    if n_taxa <= 1:
        return np.zeros(n_taxa, dtype=int)
    min_cs = max(1, int(min_cluster_size))
    n_clusters = max(1, n_taxa // min_cs)
    n_clusters = min(n_clusters, n_taxa)
    if max_cluster_size and n_clusters > 0:
        min_needed = int(np.ceil(n_taxa / max(int(max_cluster_size), 1)))
        n_clusters = max(n_clusters, min_needed)
        n_clusters = min(n_clusters, n_taxa)
    if n_clusters <= 1:
        return np.zeros(n_taxa, dtype=int)
    corr = _correlation_matrix(np.log1p(np.maximum(values, 0.0)))
    dist = np.clip(1.0 - np.abs(corr), 0.0, 2.0)
    np.fill_diagonal(dist, 0.0)
    from sklearn.cluster import AgglomerativeClustering

    model = AgglomerativeClustering(
        n_clusters=int(n_clusters), metric="precomputed", linkage="average"
    )
    try:
        return np.asarray(model.fit_predict(dist), dtype=int)
    except Exception:
        return np.zeros(n_taxa, dtype=int)


def _pairwise_glm(values: np.ndarray, min_value: float = 0.0):
    n = values.shape[0]
    fits = [[("zero",) for _ in range(n)] for _ in range(n)]
    rsq = np.zeros((n, n), dtype=float)
    prob = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                rsq[i, j] = 1.0
                prob[i, j] = 1.0
                continue
            fit = _fit_gaussian_glm(values[i], values[j], min_value)
            fits[i][j] = fit
            rsq[i, j] = _glm_rsq(fit)
            prob[i, j] = _oriented_prob(values[i], values[j], min_value)
    return fits, rsq, prob


def _which_max_from_to(rsq: np.ndarray, done: List[int]) -> Optional[Tuple[int, int]]:
    """R ``which.max.coord``: strongest edge from generated nodes to remaining."""
    n = rsq.shape[0]
    done_set = set(done)
    best = None
    best_val = 0.0
    for i in done:
        for j in range(n):
            if j in done_set:
                continue
            val = float(rsq[i, j])
            if val > best_val:
                best_val = val
                best = (i, j)
    return best


def _generate_cluster(
    values: np.ndarray,
    members: List[int],
    init_local: int,
    init_level: float,
    rng: np.random.Generator,
    noise_scale: float,
) -> Dict[int, float]:
    sub = values[members, :]
    fits, rsq, prob = _pairwise_glm(sub)
    n = len(members)
    abundance = np.zeros(n, dtype=float)
    presence = np.zeros(n, dtype=float)
    abundance[init_local] = max(0.0, float(init_level))
    presence[init_local] = 1.0
    if n == 1:
        return {members[0]: float(abundance[0])}
    edf = prob[prob > 0]
    if edf.size == 0:
        edf = np.array([0.0])
    remaining = [j for j in range(n) if j != init_local]
    draws = rng.choice(edf, size=len(remaining), replace=True)
    for j, thresh in zip(remaining, draws):
        presence[j] = 1.0 if prob[init_local, j] > float(thresh) else 0.0
    done = [init_local] + [j for j in remaining if presence[j] == 0.0]
    rsq_walk = rsq.copy()
    for j in range(n):
        if presence[j] == 0.0:
            rsq_walk[j, :] = 0.0
            rsq_walk[:, j] = 0.0
    while len(done) < n:
        edge = _which_max_from_to(rsq_walk, done)
        if edge is None:
            break
        src, dst = edge
        yval = _predict_gaussian_glm(fits[src][dst], abundance[src])
        yval = max(0.0, yval * (1.0 + float(rng.normal(0.0, noise_scale))))
        if dst in done:
            break
        abundance[dst] = yval
        done.append(dst)
    return {members[i]: float(abundance[i]) for i in range(n)}


def _glm_boil_one(
    values: np.ndarray,
    rng: np.random.Generator,
    min_cluster_size: int = 2,
    max_cluster_size: int = 100,
    noise_scale: float = 0.15,
) -> np.ndarray:
    """One synthetic sample following the R ``samovar_boil`` cluster walk."""
    n_taxa, n_samples = values.shape
    out = np.zeros(n_taxa, dtype=float)
    if n_taxa == 0 or n_samples == 0:
        return out
    labels = _cluster_taxa(values, min_cluster_size, max_cluster_size)
    clusters: Dict[int, List[int]] = {}
    for idx, lab in enumerate(labels):
        clusters.setdefault(int(lab), []).append(idx)
    cluster_ids = sorted(clusters)
    summaries = []
    for cid in cluster_ids:
        members = clusters[cid]
        summaries.append(values[members, :].mean(axis=0) if members else np.zeros(n_samples))
    summary_mat = np.vstack(summaries) if summaries else np.zeros((0, n_samples))
    _, inter_rsq, inter_prob = _pairwise_glm(summary_mat) if len(cluster_ids) else (None, None, None)

    present = np.where(values.sum(axis=1) > 0)[0]
    if present.size == 0:
        return out
    init_sp = int(rng.choice(present))
    pos = values[init_sp, values[init_sp] > 0]
    init_level = float(rng.choice(pos)) if pos.size else 0.0
    init_cluster = int(labels[init_sp])

    cluster_todo = {cid: 0 for cid in cluster_ids}
    cluster_todo[init_cluster] = 1
    if inter_prob is not None and len(cluster_ids) > 1:
        init_ci = cluster_ids.index(init_cluster)
        edf = inter_prob[inter_prob > 0]
        if edf.size == 0:
            edf = np.array([0.0])
        others = [i for i, cid in enumerate(cluster_ids) if cid != init_cluster]
        draws = rng.choice(edf, size=len(others), replace=True)
        for ci, thresh in zip(others, draws):
            cluster_todo[cluster_ids[ci]] = int(inter_prob[init_ci, ci] > float(thresh))

    generated: set = {cid for cid, flag in cluster_todo.items() if flag == 0}
    current_init = init_sp
    current_level = init_level
    current_cluster = init_cluster

    while True:
        todo = [c for c in cluster_ids if c not in generated]
        if not todo:
            break
        if current_cluster not in todo:
            current_cluster = todo[0]
            members = clusters[current_cluster]
            current_init = members[0]
            current_level = float(np.mean(values[current_init]))
        members = clusters[current_cluster]
        init_local = members.index(current_init) if current_init in members else 0
        part = _generate_cluster(
            values, members, init_local, current_level, rng, noise_scale
        )
        for sp, val in part.items():
            out[sp] = val
        generated.add(current_cluster)
        todo = [c for c in cluster_ids if c not in generated]
        if not todo or inter_rsq is None:
            break
        done_idx = [cluster_ids.index(c) for c in cluster_ids if c in generated]
        rsq_walk = inter_rsq.copy()
        for ci, cid in enumerate(cluster_ids):
            if cluster_todo.get(cid, 1) == 0:
                rsq_walk[ci, :] = 0.0
                rsq_walk[:, ci] = 0.0
        edge = _which_max_from_to(rsq_walk, done_idx)
        if edge is None:
            current_cluster = todo[0]
            members = clusters[current_cluster]
            current_init = members[0]
            current_level = float(np.mean(values[current_init]))
            continue
        src_i, dst_i = edge
        src_cid = cluster_ids[src_i]
        dst_cid = cluster_ids[dst_i]
        if dst_cid in generated:
            current_cluster = todo[0]
            continue
        src_mean = float(np.mean([out[s] for s in clusters[src_cid]]))
        src_sum = summary_mat[src_i]
        dst_sum = summary_mat[dst_i]
        fit = _fit_gaussian_glm(src_sum, dst_sum)
        cluster_level = _predict_gaussian_glm(fit, src_mean)
        dst_members = clusters[dst_cid]
        dst_sub = values[dst_members, :]
        col_ok = np.where(dst_sub.sum(axis=0) > 0)[0]
        if col_ok.size:
            pick_col = int(rng.choice(col_ok))
            row_ok = np.where(dst_sub[:, pick_col] > 0)[0]
            local = int(rng.choice(row_ok)) if row_ok.size else 0
        else:
            local = 0
        means = dst_sub.mean(axis=1)
        mean_all = float(np.mean(means)) if means.size else 1.0
        current_init = dst_members[local]
        current_level = (
            float(means[local] / mean_all * cluster_level) if mean_all else cluster_level
        )
        current_cluster = dst_cid
    return np.maximum(out, 0.0)
