"""Regenerate metagenome abundance tables from annotations or abundance CSVs.

Minimal contract: ``taxid`` × ``N_<sample>`` in, same shape out. Long
annotation tables (``taxID_*``) and phyloseq-style OTU tables are converted
first. Sequence regenerate (ISS) consumes that same abundance shape.

``config["max_genomes"]`` (CLI ``--max-genomes``, default ``inf``) caps how
many taxa are kept (top by abundance) before table regenerate and how many
non-host genomes a metagenome generator uses.

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

import math
import os
import re
import zlib
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

UNLIMITED_MAX_GENOMES = float("inf")


def parse_max_genomes(
    value: Any = None,
    *,
    default_from_env: bool = True,
) -> float:
    """Count cap for taxa / non-host genomes. ``0`` / ``inf`` / empty → unlimited."""
    if value is None and default_from_env:
        value = os.environ.get("SAMOVAR_MAX_GENOMES", "")
    if value is None or value is False:
        return UNLIMITED_MAX_GENOMES
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"", "inf", "infinity", "+inf", "none", "unlimited"}:
            return UNLIMITED_MAX_GENOMES
        try:
            value = float(text)
        except ValueError:
            return UNLIMITED_MAX_GENOMES
    try:
        number = float(value)
    except (TypeError, ValueError):
        return UNLIMITED_MAX_GENOMES
    if number <= 0 or math.isinf(number):
        return UNLIMITED_MAX_GENOMES
    return float(int(number))


def argparse_max_genomes(text: str) -> float:
    return parse_max_genomes(text, default_from_env=False)


def finite_max_genomes(
    value: Any = None,
    *,
    default_from_env: bool = False,
) -> Optional[int]:
    """``None`` when unlimited, else a positive integer cap."""
    parsed = parse_max_genomes(value, default_from_env=default_from_env)
    if parsed == UNLIMITED_MAX_GENOMES:
        return None
    return int(parsed)


def max_genomes_from_config(config: Optional[Dict[str, Any]] = None) -> float:
    cfg = config or {}
    if "max_genomes" not in cfg:
        return parse_max_genomes(None, default_from_env=True)
    return parse_max_genomes(cfg.get("max_genomes"), default_from_env=False)


def cap_sequence(items: Sequence[Any], max_genomes: Any = None) -> List[Any]:
    """Keep the first ``max_genomes`` items; ``inf`` leaves the sequence unchanged."""
    limit = finite_max_genomes(max_genomes, default_from_env=False)
    seq = list(items)
    if limit is None:
        return seq
    return seq[:limit]


def cap_matrix_taxa(matrix: pd.DataFrame, max_genomes: Any = None) -> pd.DataFrame:
    """Keep the top ``max_genomes`` rows by row-sum (taxid × sample matrix)."""
    limit = finite_max_genomes(max_genomes, default_from_env=False)
    if limit is None or matrix is None or matrix.empty or len(matrix) <= limit:
        return matrix
    totals = matrix.sum(axis=1)
    keep = totals.sort_values(ascending=False, kind="mergesort").index[:limit]
    return matrix.loc[keep]


def cap_abundance_table(frame: pd.DataFrame, max_genomes: Any = None) -> pd.DataFrame:
    """Keep the top ``max_genomes`` taxa of a ``taxid`` + ``N_*`` table."""
    from samovar.abundance import n_sample_columns

    limit = finite_max_genomes(max_genomes, default_from_env=False)
    if limit is None or frame is None or frame.empty:
        return frame
    if "taxid" not in frame.columns:
        return cap_matrix_taxa(frame, max_genomes)
    n_cols = n_sample_columns(frame)
    if not n_cols:
        return frame.iloc[:limit].copy()
    totals = frame[n_cols].sum(axis=1)
    keep = totals.sort_values(ascending=False, kind="mergesort").index[:limit]
    return frame.loc[keep].reset_index(drop=True)


def cap_abundance_tables(
    tables: Dict[str, pd.DataFrame],
    max_genomes: Any = None,
) -> Dict[str, pd.DataFrame]:
    return {name: cap_abundance_table(table, max_genomes) for name, table in tables.items()}


def cap_generate_genome_rows(
    rows: Sequence[Dict[str, str]],
    max_genomes: Any = None,
) -> List[Dict[str, str]]:
    """Keep host rows; cap non-host genomes to ``max_genomes``."""
    limit = finite_max_genomes(max_genomes, default_from_env=False)
    listed = list(rows)
    if limit is None:
        return listed
    hosts = [row for row in listed if str(row.get("host") or "") == "1"]
    meta = [row for row in listed if str(row.get("host") or "") != "1"]
    return meta[:limit] + hosts


DIRECT_MODES = frozenset(
    {"direct", "preserve", "none", "exact", "raw", "off", "false", ""}
)
GENERATIVE_MODES = frozenset({"glm", "bootstrap", "vae"})
SAMOVAR_R_MODES = frozenset({"samovar", "r", "boil"})
CAMISIM_TABLE_MODES = frozenset(
    {"camisim", "camisim-table", "camisim_table", "cami"}
)
SPARSEDOSSA2_MODE_ALIASES = {
    "sparsedossa2": "sparsedossa2-fit",
    "sparsedossa2-fit": "sparsedossa2-fit",
    "sparsedossa2_fit": "sparsedossa2-fit",
    "sd2": "sparsedossa2-fit",
    "sd2-fit": "sparsedossa2-fit",
    "sparsedossa2-stool": "sparsedossa2-stool",
    "sparsedossa2_stool": "sparsedossa2-stool",
    "sd2-stool": "sparsedossa2-stool",
    "sparsedossa2-vaginal": "sparsedossa2-vaginal",
    "sparsedossa2_vaginal": "sparsedossa2-vaginal",
    "sd2-vaginal": "sparsedossa2-vaginal",
    "sparsedossa2-ibd": "sparsedossa2-ibd",
    "sparsedossa2_ibd": "sparsedossa2-ibd",
    "sd2-ibd": "sparsedossa2-ibd",
}


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
    sd2 = SPARSEDOSSA2_MODE_ALIASES.get(low.replace("_", "-")) or SPARSEDOSSA2_MODE_ALIASES.get(low)
    if sd2:
        return "builtin", sd2
    return "custom", key


def normalize_regeneration_mode(mode: Optional[str]) -> str:
    """Built-in modes only. Custom names must go through ``resolve_regeneration_mode``."""
    kind, name = resolve_regeneration_mode(mode)
    if kind == "builtin":
        return name
    raise ValueError(
        f"Unknown regeneration_mode={mode!r}. "
        "Use direct, bootstrap, vae, glm, camisim-table, sparsedossa2-fit, "
        "sparsedossa2-stool, samovar, or an imported "
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
    from samovar.abundance import taxid_annotation_columns

    return taxid_annotation_columns(df)


def _annotator_name(col: str) -> str:
    from samovar.abundance import annotator_name

    return annotator_name(col)


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
    max_genomes: Any = None,
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
        return cap_matrix_taxa(mat, max_genomes)
    if threshold_amount and threshold_amount > 0:
        totals = mat.sum(axis=0).replace(0, 1)
        fracs = mat / totals
        keep = fracs.max(axis=1) >= float(threshold_amount)
        if not bool(keep.any()):
            keep = mat.max(axis=1) > 0
        mat = mat.loc[keep]
    return cap_matrix_taxa(mat, max_genomes)


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


def _source_matrices(data: Any) -> List[Tuple[str, pd.DataFrame]]:
    from samovar.abundance import abundance_to_matrix, input_to_abundance_tables

    out: List[Tuple[str, pd.DataFrame]] = []
    for name, table in input_to_abundance_tables(data).items():
        mat = abundance_to_matrix(table)
        if mat is not None and not mat.empty:
            out.append((str(name), mat))
    return out


def regenerate_preserve(
    data: pd.DataFrame,
    n_reads: Optional[int] = None,
    rescale: bool = False,
    threshold_amount: float = 0.0,
    max_genomes: Any = None,
) -> Dict[str, pd.DataFrame]:
    """Keep observed abundances (annotation counts or an abundance table)."""
    result: Dict[str, pd.DataFrame] = {}
    for name, mat in _source_matrices(data):
        mat = _filter_taxa(
            mat, threshold_amount, n_reads, rescale=rescale, max_genomes=max_genomes
        )
        mat = _apply_abundance_scale(mat, n_reads, rescale)
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
    max_genomes: Any = None,
) -> Dict[str, pd.DataFrame]:
    """Noisy copies of observed samples (Dirichlet-multinomial resampling).

    Each synthetic column is a resample of a real sample with sampling error,
    not a mixture of all columns. Named samples prefer the matching source.
    A fixed seed is reproducible but does not clone every column: each sample
    name gets its own RNG stream.
    """
    seed = coerce_seed(seed)
    result: Dict[str, pd.DataFrame] = {}
    for annotator, mat in _source_matrices(data):
        mat = cap_matrix_taxa(mat, max_genomes)
        if mat.shape[1] == 0:
            continue
        n_src = mat.shape[1]
        src_names = [str(c) for c in mat.columns]
        src_index = {name: i for i, name in enumerate(src_names)}
        names = synthetic_sample_names(src_names, n_samples)
        values = mat.to_numpy(dtype=float)
        synth = {}
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
        synth_mat = _filter_taxa(
            synth_mat,
            threshold_amount,
            n_reads,
            rescale=rescale,
            max_genomes=max_genomes,
        )
        synth_mat = _apply_abundance_scale(synth_mat, n_reads, rescale)
        result[annotator] = _abundance_table_from_matrix(synth_mat, annotator)
    return result


def regenerate_vae(
    data: pd.DataFrame,
    n_samples: int,
    n_reads: Optional[int] = 1000,
    threshold_amount: float = 1e-5,
    latent_dim: int = 4,
    seed: Optional[int] = 42,
    rescale: bool = False,
    max_genomes: Any = None,
) -> Dict[str, pd.DataFrame]:
    """Sample synthetic profiles with a latent linear generative model (FA)."""
    from sklearn.decomposition import FactorAnalysis

    seed = coerce_seed(seed)
    result: Dict[str, pd.DataFrame] = {}
    for annotator, mat in _source_matrices(data):
        mat = _filter_taxa(
            mat, threshold_amount, n_reads, rescale=rescale, max_genomes=max_genomes
        )
        names = synthetic_sample_names(list(mat.columns), n_samples)
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
    max_genomes: Any = None,
) -> Dict[str, pd.DataFrame]:
    """Cluster + Gaussian GLM walk matching ``samovar_boil`` / ``concotion_pour``.

    Taxa are clustered, ``y ~ x`` Gaussian GLMs are fit on co-occurring samples
    (oriented co-occurrence, Rsq edge weights), then each synthetic community
    is generated by walking that graph from a random init taxon.
    """
    seed = coerce_seed(seed)
    result: Dict[str, pd.DataFrame] = {}
    for annotator, mat in _source_matrices(data):
        mat = _filter_taxa(
            mat, threshold_amount, n_reads, rescale=rescale, max_genomes=max_genomes
        )
        if mat.empty:
            continue
        names = synthetic_sample_names(list(mat.columns), n_samples)
        values = np.maximum(mat.to_numpy(dtype=float), 0.0)
        synth_cols = {}
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
    max_genomes: Any = None,
) -> Dict[str, pd.DataFrame]:
    """CAMISIM-style log-normal abundances for taxIDs seen in the annotation table."""
    from samovar.camisim import design_communities

    seed = coerce_seed(seed)
    depth = int(n_reads) if n_reads not in (None, "", 0) else 1000
    result: Dict[str, pd.DataFrame] = {}
    for annotator, mat in _source_matrices(data):
        mat = _filter_taxa(
            mat, threshold_amount, n_reads, rescale=False, max_genomes=max_genomes
        )
        if mat.empty:
            continue
        taxids = [str(t) for t in mat.index if str(t) not in {"0", "nan", "None", ""}]
        if not taxids:
            continue
        names = synthetic_sample_names(list(mat.columns), n_samples)
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
    from samovar.abundance import load_table_input
    from samovar.table_regenerators import (
        attach_regenerator_flags,
        get_table_regenerator,
        load_samples_metadata,
    )

    source = load_table_input(annotation_dir, data)
    mode = cfg.get("table_reads_generator") or cfg.get("regeneration_mode", "direct")
    local = dict(cfg)
    local["annotation_dir"] = str(annotation_dir)
    local["output_dir"] = str(output_dir)
    local = attach_regenerator_flags(mode, local)
    regenerator = get_table_regenerator(mode)
    tables = regenerator.run(
        source,
        load_samples_metadata(local),
        local,
    )
    tables = cap_abundance_tables(tables, max_genomes_from_config(local))
    _write_abundance_tables(output_dir, tables)
    return tables


def regenerate_annotation_tables(
    annotation_dir: Union[str, Path],
    output_dir: Union[str, Path],
    config: Optional[Dict[str, Any]] = None,
    data: Optional[pd.DataFrame] = None,
    select_best: bool = True,
) -> Dict[str, pd.DataFrame]:
    """Regenerate per-annotator abundance CSVs from an annotation directory.

    If ``table_reads_generators`` names more than one method and
    ``select_best`` is true, each candidate is scored and only the winner is
    written to ``output_dir``. With ``select_best=False`` every method is kept
    under ``.table_candidates/`` and the first success is copied to ``output_dir``.
    """
    from samovar.table_regenerators import canonical_regeneration_modes
    from samovar.table_scorers import (
        rank_methods_per_annotator,
        write_table_score_plots,
        write_table_selection,
        table_score_plot_dirs,
    )

    cfg = dict(config or {})
    if data is None:
        from samovar.abundance import load_table_input

        data = load_table_input(annotation_dir)
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
    tables_by_mode: Dict[str, Dict[str, pd.DataFrame]] = {}
    generate_errors: Dict[str, str] = {}
    for mode in modes:
        cand = candidates_root / re.sub(r"[^A-Za-z0-9._-]+", "_", mode)
        one = dict(cfg)
        one["table_reads_generator"] = mode
        one["regeneration_mode"] = mode
        one["table_reads_generators"] = [mode]
        one["output_dir"] = str(cand)
        try:
            tables_by_mode[mode] = _regenerate_one_mode(annotation_dir, cand, one, data)
        except Exception as exc:
            tables_by_mode[mode] = {}
            generate_errors[mode] = str(exc)

    mixed: Dict[str, pd.DataFrame] = {}
    for mode in modes:
        if tables_by_mode.get(mode):
            mixed = tables_by_mode[mode]
            break
    if not mixed:
        raise RuntimeError(
            "No table_reads_generator produced abundance tables. "
            f"Tried: {', '.join(modes)}"
            + (f" Errors: {generate_errors}" if generate_errors else "")
        )

    if not select_best:
        _write_abundance_tables(out, mixed)
        return mixed

    from samovar.table_scorers import rank_methods_per_annotator

    try:
        ranked = rank_methods_per_annotator(
            data, tables_by_mode, scorer, config=cfg, modes=modes
        )
    except Exception as exc:
        ranked = {
            "scorer": str(scorer),
            "winner": modes[0],
            "winner_by_annotator": {},
            "by_annotator": {},
            "candidates": [
                {
                    "ok": False,
                    "error": str(exc),
                    "rank_value": 1.0e9,
                    "scorer": str(scorer),
                }
            ],
            "tables": {},
        }
    for mode, err in generate_errors.items():
        ranked.setdefault("generate_errors", {})[mode] = err
    mixed = ranked.get("tables") or {}
    if not mixed:
        for mode in modes:
            if tables_by_mode.get(mode):
                mixed = tables_by_mode[mode]
                ranked["winner"] = mode
                ranked["tables"] = mixed
                break
    if not mixed:
        raise RuntimeError(
            "No table_reads_generator produced abundance tables. "
            f"Tried: {', '.join(modes)}"
            + (f" Errors: {generate_errors}" if generate_errors else "")
        )
    _write_abundance_tables(out, mixed)
    payload = {k: v for k, v in ranked.items() if k != "tables"}
    write_table_selection(out / "table_selection.json", payload)
    write_table_selection(candidates_root / "table_selection.json", payload)
    try:
        for dest in table_score_plot_dirs(out):
            write_table_score_plots(payload, dest, config=cfg)
    except Exception:
        pass
    return mixed


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
        "max_genomes": UNLIMITED_MAX_GENOMES,
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
    from samovar.abundance import observed_n_samples

    observed = max(int(observed_n_samples(data)), 1)
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


def stage_regenerate_tables(
    output_dir: Union[str, Path],
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, pd.DataFrame]:
    """Checkpoint ``regenerate_tables``: observed abundance → regenerated CSVs."""
    from samovar.abundance import (
        has_abundance_tables,
        materialize_observed_abundance,
        observed_abundance_dir,
        regenerated_abundance_dir,
    )

    root = Path(output_dir)
    observed = observed_abundance_dir(root)
    if not has_abundance_tables(observed):
        materialize_observed_abundance(root)
    dest = regenerated_abundance_dir(root)
    return regenerate_annotation_tables(
        observed, dest, dict(config or {}), select_best=False
    )


def main(argv: Optional[List[str]] = None) -> int:
    import argparse

    import yaml

    parser = argparse.ArgumentParser(prog="python -m samovar.regenerate")
    sub = parser.add_subparsers(dest="command", required=True)
    tables = sub.add_parser("tables", help="Regenerate abundance tables for a run directory")
    tables.add_argument("--output_dir", "--outdir", dest="output_dir", required=True)
    tables.add_argument("--config", default="", help="YAML (annotation2iss / table regen)")
    args = parser.parse_args(argv)
    if args.command == "tables":
        cfg: Dict[str, Any] = {}
        if args.config:
            cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
        written = stage_regenerate_tables(args.output_dir, cfg)
        print(f"Wrote {len(written)} regenerated abundance table(s)")
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

