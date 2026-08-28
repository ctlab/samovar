"""Score regenerated abundance tables against the observed annotation.

Built-ins at this stage (``--table-score`` / ``table_score``) share the same
ranking block: one row per regeneration method, lower ``rank_value`` wins.

- ``shannon_ks``: KS on per-sample Shannon entropy
- ``bray_ks``: KS on pairwise Bray–Curtis distances
- ``sparsedossa2_cv``: SparseDOSSA2 ``fitCV`` (higher GOF → lower rank_value)

Imported ``--type table-scoring`` Python modules fill the same ranking block
via ``score_annotator(...)`` or ``score_table(observed, generated, config)``.
Annotation/viz ``--type scoring`` (``score(inputs, dest, config)``) is a
different contract and is not used here.
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
        "bray_ks",
        "bray",
        "bray_curtis_ks",
        "sparsedossa2_cv",
        "sparsedossa2-cv",
        "sparsedossa2",
        "sd2_cv",
        "sd2-cv",
        "fitcv",
    }
)

PathLike = Union[str, Path]
TABLE_SCORING_GROUP = "table_scoring"


class MissingTableScorerError(ValueError):
    """``tools.<name>`` is missing or is not ``--type table-scoring``."""


def flags_apply_to_table_scorer(target: str, name: Optional[str] = None) -> bool:
    from samovar.main_config import flags_target_matches

    names = [name] if name else []
    return flags_target_matches(
        target,
        *names,
        groups=("table_scoring", "table-scoring", "table_score"),
    )


def lookup_table_scorer(name: str):
    """Return ``(config_key, parsed_spec)`` for an imported table scorer."""
    from samovar.main_config import iter_tools, parse_tool_entry, tool_path
    from samovar.paths import load_config

    key = str(name or "").strip()
    if not key:
        raise MissingTableScorerError(
            "Empty table scorer name. Import a tool with "
            "`samovar tools import -n NAME --type table-scoring`."
        )
    tools = iter_tools(load_config())
    spec = tools.get(key)
    matched = key
    if spec is None:
        low = key.lower()
        for stored, value in tools.items():
            if stored.lower() == low:
                spec = value
                matched = stored
                break
    if spec is None:
        raise MissingTableScorerError(
            f"table-scoring tool {key!r} is not in the install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type table-scoring`."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != TABLE_SCORING_GROUP:
        raise MissingTableScorerError(
            f"tools.{matched} has group {group!r}, expected {TABLE_SCORING_GROUP!r}. "
            "Re-import with --type table-scoring (not --type scoring)."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingTableScorerError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return matched, parsed


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
    """Shannon diversity per (annotator, sample) from annotation or abundance."""
    from samovar.abundance import input_to_abundance_tables

    chunks = [
        shannon_vector_from_table(table)
        for table in input_to_abundance_tables(data).values()
    ]
    nonempty = [c for c in chunks if c.size]
    if not nonempty:
        return np.array([], dtype=float)
    return np.concatenate(nonempty)


def shannon_vector_from_table(table: pd.DataFrame) -> np.ndarray:
    from samovar.abundance import abundance_to_matrix, is_abundance_table, normalize_abundance_table

    if table is None or getattr(table, "empty", True):
        return np.array([], dtype=float)
    if is_abundance_table(table) or "taxid" in getattr(table, "columns", []):
        work = normalize_abundance_table(table)
        mat = _drop_unclassified_rows(abundance_to_matrix(work))
    else:
        return np.array([], dtype=float)
    if mat.empty:
        return np.array([], dtype=float)
    values: List[float] = []
    for col in mat.columns:
        values.append(shannon_entropy(mat[col].to_numpy(dtype=float)))
    return np.asarray(values, dtype=float)


def shannon_vector_from_tables(tables: Dict[str, pd.DataFrame]) -> np.ndarray:
    """Shannon diversity per (annotator table, N_* sample column)."""
    chunks = [shannon_vector_from_table(table) for table in (tables or {}).values()]
    nonempty = [c for c in chunks if c.size]
    if not nonempty:
        return np.array([], dtype=float)
    return np.concatenate(nonempty)


def shannon_vector_from_annotation_annotator(data: pd.DataFrame, annotator: str) -> np.ndarray:
    from samovar.abundance import input_to_abundance_tables

    tables = input_to_abundance_tables(data)
    if annotator in tables:
        return shannon_vector_from_table(tables[annotator])
    if len(tables) == 1:
        return shannon_vector_from_table(next(iter(tables.values())))
    return np.array([], dtype=float)


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
    key = str(name or "shannon_ks").strip()
    low = key.lower().replace("-", "_")
    if low in {"sparsedossa2_cv", "sparsedossa2", "sd2_cv", "fitcv"}:
        return "sparsedossa2_cv"
    if low in {"shannon_ks", "ks", "shannon", "alpha_ks"}:
        return "shannon_ks"
    if low in {"bray_ks", "bray", "bray_curtis_ks"}:
        return "bray_ks"
    if not key:
        return "shannon_ks"
    try:
        matched, _spec = lookup_table_scorer(key)
        return matched
    except MissingTableScorerError:
        pass
    raise ValueError(
        f"Unknown table scorer {name!r}. Built-in: shannon_ks, sparsedossa2_cv, bray_ks. "
        "Or import one with `samovar tools import --type table-scoring` and pass its name to --table-score."
    )


def _tables_to_matrices(tables: Dict[str, pd.DataFrame]) -> Dict[str, pd.DataFrame]:
    from samovar.abundance import abundance_to_matrix, is_abundance_table, normalize_abundance_table

    matrices: Dict[str, pd.DataFrame] = {}
    for name, table in (tables or {}).items():
        if table is None or getattr(table, "empty", True):
            continue
        if is_abundance_table(table) or "taxid" in table.columns:
            mat = abundance_to_matrix(normalize_abundance_table(table))
        else:
            continue
        if mat is None or mat.empty:
            continue
        matrices[str(name)] = mat
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


def _bray_vector_from_tables(tables: Dict[str, pd.DataFrame]) -> np.ndarray:
    chunks = [pairwise_bray_curtis(table) for table in (tables or {}).values()]
    nonempty = [c for c in chunks if c.size]
    if not nonempty:
        return np.array([], dtype=float)
    return np.concatenate(nonempty)


def pairwise_bray_curtis(table: pd.DataFrame) -> np.ndarray:
    """Upper-triangle Bray–Curtis distances between sample columns."""
    from samovar.abundance import abundance_to_matrix, normalize_abundance_table

    if table is None or getattr(table, "empty", True):
        return np.array([], dtype=float)
    mat = _drop_unclassified_rows(abundance_to_matrix(normalize_abundance_table(table)))
    if mat.shape[1] < 2:
        return np.array([], dtype=float)
    arr = np.maximum(mat.to_numpy(dtype=float), 0.0)
    dists: List[float] = []
    n = arr.shape[1]
    for i in range(n):
        for j in range(i + 1, n):
            num = float(np.abs(arr[:, i] - arr[:, j]).sum())
            den = float(arr[:, i].sum() + arr[:, j].sum())
            dists.append(num / den if den else 0.0)
    return np.asarray(dists, dtype=float)


def ks_bray(observed: pd.DataFrame, generated: pd.DataFrame) -> Dict[str, Any]:
    """KS test on pairwise Bray–Curtis distances (lower D is a closer community)."""
    row = ks_shannon(pairwise_bray_curtis(observed), pairwise_bray_curtis(generated))
    row["scorer"] = "bray_ks"
    return row


def score_bray_ks_annotator(
    annotation: Any,
    tables_by_mode: Dict[str, Dict[str, pd.DataFrame]],
    annotator: str,
    modes: Sequence[str],
) -> Dict[str, Any]:
    from samovar.abundance import input_to_abundance_tables

    observed_tables = input_to_abundance_tables(annotation)
    observed = observed_tables.get(annotator)
    if observed is None and len(observed_tables) == 1:
        observed = next(iter(observed_tables.values()))
    generated_by_mode: Dict[str, pd.DataFrame] = {}
    rows: List[Dict[str, Any]] = []
    for mode in modes:
        table = (tables_by_mode.get(mode) or {}).get(annotator)
        generated_by_mode[mode] = table if table is not None else pd.DataFrame()
        score = ks_bray(observed if observed is not None else pd.DataFrame(), generated_by_mode[mode])
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
                pair = ks_bray(generated_by_mode.get(mode_i), generated_by_mode.get(mode_j))
                grid[mode_i][mode_j] = pair.get("ks_statistic")
    return {
        "scorer": "bray_ks",
        "annotator": annotator,
        "winner": pick_best_table_method(rows),
        "score_matrix": grid,
        "cv_matrix": grid,
        "candidates": rows,
        "metric_name": "ks_statistic",
    }


def _load_imported_table_scorer(name: str) -> Any:
    import importlib.util

    from samovar.main_config import tool_path

    matched, spec = lookup_table_scorer(name)
    path = Path(tool_path(spec, matched)).expanduser()
    if path.suffix.lower() != ".py" or not path.is_file():
        raise ValueError(
            f"Imported table scorer {name!r} must be a Python module with "
            "score_annotator() or score_table() "
            "(annotation --type scoring is a different contract)."
        )
    loaded = importlib.util.spec_from_file_location(f"samovar_table_scorer_{matched}", path)
    if loaded is None or loaded.loader is None:
        raise ValueError(f"Cannot load table scorer {path}")
    module = importlib.util.module_from_spec(loaded)
    loaded.loader.exec_module(module)
    return module


def score_imported_annotator(
    name: str,
    annotation: Any,
    tables_by_mode: Dict[str, Dict[str, pd.DataFrame]],
    annotator: str,
    modes: Sequence[str],
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Same ranking block as shannon_ks, filled by an imported --type table-scoring module."""
    module = _load_imported_table_scorer(name)
    cfg = dict(config or {})
    fn_ann = getattr(module, "score_annotator", None)
    if callable(fn_ann):
        block = fn_ann(annotation, tables_by_mode, annotator, list(modes), cfg)
        if not isinstance(block, dict):
            raise TypeError(f"{name}.score_annotator() must return a dict")
        block.setdefault("scorer", name)
        block.setdefault("annotator", annotator)
        block.setdefault("metric_name", "rank_value")
        return block
    fn_one = getattr(module, "score_table", None)
    if not callable(fn_one):
        raise ValueError(
            f"Imported scorer {name!r} needs score_annotator(...) or score_table(observed, generated, config). "
            "The annotation-stage score(inputs, dest, config) hook is not used for table regenerate."
        )
    from samovar.abundance import input_to_abundance_tables

    observed_tables = input_to_abundance_tables(annotation)
    observed = observed_tables.get(annotator)
    if observed is None and len(observed_tables) == 1:
        observed = next(iter(observed_tables.values()))
    rows: List[Dict[str, Any]] = []
    for mode in modes:
        generated = (tables_by_mode.get(mode) or {}).get(annotator)
        raw = fn_one(
            observed if observed is not None else pd.DataFrame(),
            generated if generated is not None else pd.DataFrame(),
            cfg,
        )
        if not isinstance(raw, dict):
            raise TypeError(f"{name}.score_table() must return a dict")
        row = dict(raw)
        row.setdefault("scorer", name)
        row["mode"] = mode
        row["annotator"] = annotator
        row.setdefault("ok", generated is not None and not getattr(generated, "empty", True))
        if "rank_value" not in row:
            row["rank_value"] = float(row.get("metric") or row.get("ks_statistic") or 1.0e9)
        rows.append(row)
    return {
        "scorer": name,
        "annotator": annotator,
        "winner": pick_best_table_method(rows),
        "score_matrix": {},
        "cv_matrix": {},
        "candidates": rows,
        "metric_name": str(rows[0].get("metric_name") or "rank_value") if rows else "rank_value",
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
        elif kind == "shannon_ks":
            block = score_shannon_ks_annotator(
                annotation, tables_by_mode, annotator, mode_list
            )
        elif kind == "bray_ks":
            block = score_bray_ks_annotator(
                annotation, tables_by_mode, annotator, mode_list
            )
        else:
            block = score_imported_annotator(
                kind, annotation, tables_by_mode, annotator, mode_list, config=config
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
    if kind == "bray_ks":
        from samovar.abundance import input_to_abundance_tables

        row = ks_shannon(
            _bray_vector_from_tables(input_to_abundance_tables(annotation)),
            _bray_vector_from_tables(tables),
        )
        row["scorer"] = "bray_ks"
        return row
    if kind == "sparsedossa2_cv":
        return score_sparsedossa2_cv(annotation, tables, config=config)
    ranked = rank_methods_per_annotator(
        annotation, {"_one": tables}, scorer=kind, config=config, modes=["_one"]
    )
    cands = list(ranked.get("candidates") or [])
    if not cands:
        raise ValueError(f"Imported table scorer {kind!r} returned no candidates")
    return min(cands, key=lambda r: float(r.get("rank_value", 1.0e9)))


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


def load_tables_by_mode_from_run(output_dir: PathLike) -> Dict[str, Dict[str, pd.DataFrame]]:
    from samovar.abundance import load_abundance_dir, regenerated_abundance_dir

    dest = regenerated_abundance_dir(output_dir)
    candidates = dest / ".table_candidates"
    tables_by_mode: Dict[str, Dict[str, pd.DataFrame]] = {}
    if candidates.is_dir():
        for child in sorted(candidates.iterdir()):
            if not child.is_dir():
                continue
            loaded = load_abundance_dir(child)
            if loaded:
                tables_by_mode[child.name] = loaded
    if not tables_by_mode:
        loaded = load_abundance_dir(dest)
        if loaded:
            tables_by_mode["generated"] = loaded
    return tables_by_mode


def stage_score_regenerated_tables(
    output_dir: PathLike,
    config: Optional[Dict[str, Any]] = None,
) -> Dict[str, Any]:
    """Checkpoint ``score_regenerated_tables``: rank methods and write winner CSVs."""
    from samovar.abundance import (
        has_abundance_tables,
        load_table_input,
        materialize_observed_abundance,
        observed_abundance_dir,
        regenerated_abundance_dir,
        write_abundance_dir,
    )

    cfg = dict(config or {})
    root = Path(output_dir)
    observed_dir = observed_abundance_dir(root)
    if not has_abundance_tables(observed_dir):
        materialize_observed_abundance(root)
    observed = load_table_input(observed_dir)
    tables_by_mode = load_tables_by_mode_from_run(root)
    if not tables_by_mode:
        raise FileNotFoundError(
            f"No regenerated abundance tables under {regenerated_abundance_dir(root)}"
        )
    scorer = cfg.get("table_score") or cfg.get("table_reads_scorer") or "shannon_ks"
    ranked = rank_methods_per_annotator(
        observed, tables_by_mode, scorer, config=cfg, modes=list(tables_by_mode)
    )
    mixed = ranked.get("tables") or {}
    dest = regenerated_abundance_dir(root)
    if mixed:
        write_abundance_dir(dest, mixed)
    payload = {k: v for k, v in ranked.items() if k != "tables"}
    write_table_selection(dest / "table_selection.json", payload)
    try:
        for plot_dest in table_score_plot_dirs(dest):
            write_table_score_plots(payload, plot_dest, config=cfg)
    except Exception:
        pass
    return ranked


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    import yaml

    parser = argparse.ArgumentParser(prog="python -m samovar.table_scorers")
    sub = parser.add_subparsers(dest="command", required=True)
    stage = sub.add_parser("stage", help="Score regenerated abundance tables for a run")
    stage.add_argument("--output_dir", "--outdir", dest="output_dir", required=True)
    stage.add_argument("--config", default="")
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.command == "stage":
        cfg: Dict[str, Any] = {}
        if args.config:
            cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
        ranked = stage_score_regenerated_tables(args.output_dir, cfg)
        print(f"scorer={ranked.get('scorer')} winner={ranked.get('winner')}")
        return 0
    return 2


if __name__ == "__main__":
    raise SystemExit(main())

