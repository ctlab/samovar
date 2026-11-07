"""Concordance of two beta-diversity distance matrices.

``mantel`` is the built-in scorer. Another scorer is a ``--type concordance``
module with ``score(left, right, config)``, or a Python path in ``concordance``.
The dispatch aligns sample ids and checks the matrices before the scorer runs.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import numpy as np
import pandas as pd

GROUP = "concordance"
_FIELDS = (
    "method",
    "version",
    "statistic",
    "pvalue",
    "n_samples",
    "seed",
    "permutations",
    "status",
    "reason",
)


class MissingConcordanceScorerError(ValueError):
    """``tools.<name>`` is missing or is not ``--type concordance``."""


def prepare_distance_pair(left, right) -> Tuple[Optional[pd.DataFrame], Optional[pd.DataFrame], Dict[str, Any]]:
    """Align on shared sample ids. Samples present on only one side are omitted."""
    checked = []
    for label, frame in (("left", left), ("right", right)):
        matrix, reason = _check_matrix(frame, label)
        if reason:
            return None, None, {"status": "invalid distance matrix", "reason": reason, "omitted_samples": []}
        checked.append(matrix)
    left_matrix, right_matrix = checked
    shared = [sample for sample in left_matrix.index if sample in set(right_matrix.index)]
    omitted = sorted(set(left_matrix.index).symmetric_difference(right_matrix.index))
    left_matrix = left_matrix.loc[shared, shared]
    right_matrix = right_matrix.loc[shared, shared]
    return left_matrix, right_matrix, {"status": None, "reason": "", "omitted_samples": omitted, "n_samples": len(shared)}


def score_concordance(left, right, config=None) -> Dict[str, Any]:
    """Score one aligned pair. ``config['concordance']`` selects the scorer."""
    cfg = dict(config or {})
    name = str(cfg.get("concordance") or "mantel").strip() or "mantel"
    left_matrix, right_matrix, prepared = prepare_distance_pair(left, right)
    base = {
        "method": name,
        "version": "",
        "statistic": None,
        "pvalue": None,
        "n_samples": prepared.get("n_samples", 0),
        "seed": cfg.get("seed", 42),
        "permutations": cfg.get("permutations", 999),
        "status": prepared["status"],
        "reason": prepared["reason"],
        "omitted_samples": prepared["omitted_samples"],
        "scorer": name,
    }
    if prepared["status"]:
        return base
    module = _load_scorer(name)
    raw = module.score(left_matrix, right_matrix, cfg)
    if not isinstance(raw, dict):
        raise TypeError(f"{name}.score() must return a dict")
    row = dict(base)
    for key in _FIELDS:
        if key in raw:
            row[key] = raw[key]
    for key, value in raw.items():
        if key not in row:
            row[key] = value
    row["scorer"] = name
    row["omitted_samples"] = prepared["omitted_samples"]
    row["n_samples"] = int(raw.get("n_samples") if raw.get("n_samples") is not None else prepared["n_samples"])
    return row


def _check_matrix(frame, label):
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return None, f"{label} distance matrix is empty"
    matrix = frame.copy()
    matrix.index = matrix.index.astype(str)
    matrix.columns = matrix.columns.astype(str)
    if set(matrix.index) != set(matrix.columns):
        return None, f"{label} sample ids differ between rows and columns"
    matrix = matrix.loc[list(matrix.index), list(matrix.index)]
    values = matrix.to_numpy(dtype=float)
    if values.shape[0] != values.shape[1] or not np.isfinite(values).all():
        return None, f"{label} distance matrix is not finite and square"
    if not np.allclose(values, values.T, atol=1e-8):
        return None, f"{label} distance matrix is not symmetric"
    if not np.allclose(np.diag(values), 0, atol=1e-8):
        return None, f"{label} distance matrix diagonal is not zero"
    return matrix, None


def _load_scorer(name: str):
    key = name.lower().replace("-", "_")
    if key == "mantel":
        from samovar.baselines import mantel
        return mantel
    path = Path(name).expanduser()
    if path.suffix.lower() == ".py" and path.is_file():
        return _load_module(path, f"samovar_concordance_{path.stem}")
    matched, spec = lookup_concordance_scorer(name)
    from samovar.main_config import tool_path
    module_path = Path(tool_path(spec, matched)).expanduser()
    if module_path.suffix.lower() != ".py" or not module_path.is_file():
        raise MissingConcordanceScorerError(
            f"Concordance scorer {name!r} must be a Python module with score(left, right, config)."
        )
    return _load_module(module_path, f"samovar_concordance_{matched}")


def lookup_concordance_scorer(name: str):
    """Return ``(config_key, parsed_spec)`` for an imported concordance scorer."""
    from samovar.main_config import iter_tools, parse_tool_entry, tool_path
    from samovar.paths import load_config

    key = str(name or "").strip()
    tools = iter_tools(load_config())
    spec = tools.get(key)
    matched = key
    if spec is None:
        for stored, value in tools.items():
            if stored.lower() == key.lower():
                spec = value
                matched = stored
                break
    if spec is None:
        raise MissingConcordanceScorerError(
            f"concordance tool {key!r} is not in the install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type concordance`."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != GROUP:
        raise MissingConcordanceScorerError(
            f"tools.{matched} has group {group!r}, expected {GROUP!r}. "
            "Re-import with --type concordance."
        )
    if not tool_path(parsed, matched):
        raise MissingConcordanceScorerError(f"tools.{matched} has an empty path.")
    return matched, parsed


def _load_module(path: Path, name: str):
    loaded = importlib.util.spec_from_file_location(name, path)
    if loaded is None or loaded.loader is None:
        raise MissingConcordanceScorerError(f"Cannot load concordance scorer {path}")
    module = importlib.util.module_from_spec(loaded)
    loaded.loader.exec_module(module)
    if not callable(getattr(module, "score", None)):
        raise MissingConcordanceScorerError(f"{path} must define score(left, right, config)")
    return module


def mantel_score(left, right, config=None):
    """Pearson correlation of the upper triangles.

    Permutations shuffle sample labels of the right matrix, not individual
    cells. The p-value is two-sided:
    ``(1 + #{|r_perm| >= |r|}) / (permutations + 1)``.
    """
    cfg = dict(config or {})
    permutations = int(cfg.get("permutations") or 999)
    seed = 42 if cfg.get("seed") is None else int(cfg["seed"])
    left_values = pd.DataFrame(left).to_numpy(dtype=float)
    right_values = pd.DataFrame(right).to_numpy(dtype=float)
    size = left_values.shape[0]
    result = {
        "method": "mantel",
        "version": "pearson-upper-triangle",
        "statistic": None,
        "pvalue": None,
        "n_samples": size,
        "seed": seed,
        "permutations": permutations,
        "status": None,
        "reason": "",
    }
    if size < 3 or right_values.shape[0] < 3:
        result["status"] = "fewer than three shared samples"
        result["reason"] = "Mantel needs at least three shared samples"
        return result
    observed_left = left_values[np.triu_indices(size, k=1)]
    observed_right = right_values[np.triu_indices(size, k=1)]
    if np.allclose(observed_left, observed_left[0]) or np.allclose(observed_right, observed_right[0]):
        result["status"] = "constant distances"
        result["reason"] = "A constant upper triangle has no Pearson correlation"
        return result
    statistic = _pearson(observed_left, observed_right)
    result["statistic"] = statistic
    generator = np.random.default_rng(seed)
    order = np.arange(size)
    hits = 0
    for _ in range(permutations):
        shuffled = generator.permutation(order)
        permuted = right_values[np.ix_(shuffled, shuffled)]
        correlation = _pearson(observed_left, permuted[np.triu_indices(size, k=1)])
        if abs(correlation) + 1e-12 >= abs(statistic):
            hits += 1
    result["pvalue"] = (hits + 1) / (permutations + 1)
    return result


def _pearson(left, right):
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    left = left - left.mean()
    right = right - right.mean()
    scale = float(np.sqrt((left * left).sum() * (right * right).sum()))
    if scale <= 0:
        return float("nan")
    return float((left * right).sum() / scale)
