"""Score regenerated abundance tables against the observed annotation.

The default builtin ``shannon_ks`` compares Shannon alpha-diversity
distributions (one value per annotator×sample) with a two-sample Kolmogorov–
Smirnov test. Lower KS statistic means the generated community better matches
the observed diversity; that candidate is selected for ISS.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Union

import numpy as np
import pandas as pd

UNCLASSIFIED = frozenset({"0", "unclassified", "unclassified_root", "nan", "none", ""})
BUILTIN_TABLE_SCORERS = frozenset({"shannon_ks", "ks", "shannon", "alpha_ks"})

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


def shannon_vector_from_tables(tables: Dict[str, pd.DataFrame]) -> np.ndarray:
    """Shannon diversity per (annotator table, N_* sample column)."""
    values: List[float] = []
    for table in (tables or {}).values():
        if table is None or table.empty or "taxid" not in table.columns:
            continue
        work = table.copy()
        work["taxid"] = work["taxid"].astype(str)
        n_cols = [c for c in work.columns if str(c).startswith("N_")]
        if not n_cols:
            continue
        work = work.loc[~work["taxid"].str.lower().isin(UNCLASSIFIED)]
        for col in n_cols:
            values.append(shannon_entropy(pd.to_numeric(work[col], errors="coerce").fillna(0)))
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
        return {
            "scorer": "shannon_ks",
            "ks_statistic": 1.0,
            "pvalue": 0.0,
            "n_observed": n_obs,
            "n_generated": n_gen,
            "rank_value": 1.0,
        }
    from scipy.stats import ks_2samp

    result = ks_2samp(obs, gen, alternative="two-sided", mode="auto")
    stat = float(result.statistic)
    pvalue = float(result.pvalue)
    if not math.isfinite(stat):
        stat = 1.0
    if not math.isfinite(pvalue):
        pvalue = 0.0
    return {
        "scorer": "shannon_ks",
        "ks_statistic": stat,
        "pvalue": pvalue,
        "n_observed": n_obs,
        "n_generated": n_gen,
        "rank_value": stat,
        "observed_shannon": [float(x) for x in obs],
        "generated_shannon": [float(x) for x in gen],
    }


def canonicalize_table_scorer(name: Optional[str]) -> str:
    key = str(name or "shannon_ks").strip().lower().replace("-", "_")
    if key in BUILTIN_TABLE_SCORERS:
        return "shannon_ks"
    if not key:
        return "shannon_ks"
    raise ValueError(
        f"Unknown table scorer {name!r}. Built-in: shannon_ks "
        "(KS test on Shannon alpha-diversity)."
    )


def score_generated_tables(
    annotation: pd.DataFrame,
    tables: Dict[str, pd.DataFrame],
    scorer: Optional[str] = None,
) -> Dict[str, Any]:
    kind = canonicalize_table_scorer(scorer)
    if kind == "shannon_ks":
        return ks_shannon(
            shannon_vector_from_annotation(annotation),
            shannon_vector_from_tables(tables),
        )
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
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path
