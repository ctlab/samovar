"""Grouped whole-row resampling of an existing regenerated CAT."""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any, Dict, Mapping, Sequence

import numpy as np
import pandas as pd

SCHEMA_VERSION = 1
SOLVER = "raking-v1"
PROVENANCE_COLUMNS = ("draw_id", "source_row", "resampling_group")
_FORBIDDEN_MATCH = {
    "true",
    "taxid_samovar",
    "seq",
    "sample",
    "draw_id",
    "source_row",
    "resampling_group",
}


def _frame(table: Any) -> pd.DataFrame:
    if isinstance(table, pd.DataFrame):
        return table.copy()
    raise TypeError("CAT resampling expects a DataFrame")


def _validate_columns(initial: pd.DataFrame, pools: Mapping[str, pd.DataFrame], columns) -> list[str]:
    selected = [str(column) for column in columns]
    if not selected:
        raise ValueError("resampling needs at least one matching column")
    if len(selected) != len(set(selected)):
        raise ValueError("resampling matching columns must be unique")
    forbidden = [column for column in selected if column.strip().lower() in _FORBIDDEN_MATCH]
    if forbidden:
        raise ValueError("resampling matching cannot use " + ", ".join(forbidden))
    missing = [column for column in selected if column not in initial.columns]
    if missing:
        raise ValueError("initial CAT is missing matching columns: " + ", ".join(missing))
    for name, pool in pools.items():
        absent = [column for column in selected if column not in pool.columns]
        if absent:
            raise ValueError(f"regenerated dataset {name!r} is missing: {', '.join(absent)}")
    return selected


def _is_taxid_column(column: str) -> bool:
    low = str(column).strip().lower()
    if "confidence" in low or low.endswith("_conf") or low.endswith("_abstain_reason"):
        return False
    return low == "taxid" or low.startswith("taxid_")


def _tokens(values: pd.Series, column: str | None = None) -> pd.Series:
    text = values.astype("string")
    if column and _is_taxid_column(column):
        numeric = pd.to_numeric(values, errors="coerce")
        integral = numeric.notna() & (np.abs(numeric % 1) < 1e-9)
        if integral.any():
            text = text.astype("object")
            text.loc[integral] = numeric.loc[integral].round().astype("int64").astype(str)
    return text.fillna("<NA>").astype(str)


def _is_categorical(column: str, target: pd.Series, pool: pd.Series) -> bool:
    low = str(column).strip().lower()
    if (
        _is_taxid_column(column)
        or low.endswith("_kmer_cluster")
        or low == "mag_id"
        or low.endswith("_mag_id")
    ):
        return True
    target_numeric = pd.to_numeric(target, errors="coerce")
    pool_numeric = pd.to_numeric(pool, errors="coerce")
    target_values = int(target.notna().sum())
    pool_values = int(pool.notna().sum())
    return (
        not target_values
        or not pool_values
        or int(target_numeric.notna().sum()) < target_values
        or int(pool_numeric.notna().sum()) < pool_values
    )


def _descriptor(column: str, target: pd.Series, pool: pd.Series) -> dict:
    categorical = _is_categorical(column, target, pool)
    if categorical:
        target_values = _tokens(target, column)
        pool_values = _tokens(pool, column)
    else:
        target_numeric = pd.to_numeric(target, errors="coerce")
        pool_numeric = pd.to_numeric(pool, errors="coerce")
        quantiles = np.linspace(0, 1, min(20, max(2, len(target_numeric))) + 1)
        edges = np.unique(target_numeric.quantile(quantiles).dropna().to_numpy(dtype=float))
        if len(edges) < 2:
            value = float(edges[0]) if len(edges) else 0.0
            epsilon = max(abs(value) * 1e-9, 1e-9)
            edges = np.asarray([-np.inf, value - epsilon, value + epsilon, np.inf])
        else:
            midpoints = (edges[:-1] + edges[1:]) / 2
            edges = np.concatenate(([-np.inf], midpoints, [np.inf]))
        target_values = pd.Series(
            pd.cut(target_numeric, edges, include_lowest=True, duplicates="drop").astype(str),
            index=target.index,
        )
        pool_values = pd.Series(
            pd.cut(pool_numeric, edges, include_lowest=True, duplicates="drop").astype(str),
            index=pool.index,
        )
    target_counts = target_values.value_counts(normalize=True)
    pool_support = set(pool_values)
    missing = [str(value) for value in target_counts.index if value not in pool_support]
    return {
        "column": column,
        "categorical": categorical,
        "target": target_values,
        "pool": pool_values,
        "target_frequency": target_counts.to_dict(),
        "missing": missing,
    }


def _joint_frame(descriptors: Sequence[dict], side: str) -> pd.DataFrame:
    return pd.concat(
        [item[side].rename(item["column"]).reset_index(drop=True) for item in descriptors],
        axis=1,
    )


def _joint_descriptor(descriptors: Sequence[dict]) -> dict:
    target_frame = _joint_frame(descriptors, "target")
    pool_frame = _joint_frame(descriptors, "pool")
    target = pd.Series(
        [int(value) for value in pd.util.hash_pandas_object(target_frame, index=False)],
        index=target_frame.index,
    )
    pool = pd.Series(
        [int(value) for value in pd.util.hash_pandas_object(pool_frame, index=False)],
        index=pool_frame.index,
    )
    labels = (
        target_frame.astype("string").fillna("<NA>").astype(str).agg(" | ".join, axis=1)
    )
    label_for_key = dict(zip(target.tolist(), labels.tolist()))
    target_counts = target.value_counts(normalize=True)
    pool_support = set(pool.tolist())
    missing = [
        label_for_key[value]
        for value in target_counts.index
        if value not in pool_support
    ]
    return {
        "column": "<joint>",
        "categorical": True,
        "target": target,
        "pool": pool,
        "target_frequency": target_counts.to_dict(),
        "missing": missing,
    }


def _rake_weights(descriptors: Sequence[dict], n_rows: int) -> np.ndarray:
    weights = np.full(n_rows, 1.0 / n_rows)
    for _ in range(25):
        previous = weights.copy()
        for descriptor in descriptors:
            pool_values = descriptor["pool"].reset_index(drop=True)
            current = pd.Series(weights).groupby(pool_values).sum()
            factors = {
                value: float(target) / float(current.get(value, 0.0))
                for value, target in descriptor["target_frequency"].items()
                if current.get(value, 0.0) > 0
            }
            adjusted = pool_values.map(factors).fillna(0.0).to_numpy(dtype=float)
            candidate = weights * adjusted
            total = float(candidate.sum())
            if total > 0:
                weights = candidate / total
        if float(np.abs(weights - previous).sum()) < 1e-10:
            break
    return weights


def _categorical_distance(left: pd.Series, right: pd.Series, column: str) -> float:
    left_freq = _tokens(left, column).value_counts(normalize=True)
    right_freq = _tokens(right, column).value_counts(normalize=True)
    keys = left_freq.index.union(right_freq.index)
    return float(
        np.abs(
            left_freq.reindex(keys, fill_value=0).to_numpy()
            - right_freq.reindex(keys, fill_value=0).to_numpy()
        ).sum()
        / 2
    )


def _numeric_distance(left: pd.Series, right: pd.Series) -> float:
    a = np.sort(pd.to_numeric(left, errors="coerce").dropna().to_numpy(dtype=float))
    b = np.sort(pd.to_numeric(right, errors="coerce").dropna().to_numpy(dtype=float))
    if not len(a) or not len(b):
        return 1.0
    values = np.unique(np.concatenate([a, b]))
    return float(
        np.max(
            np.abs(
                np.searchsorted(a, values, side="right") / len(a)
                - np.searchsorted(b, values, side="right") / len(b)
            )
        )
    )


def distribution_mismatch(left: pd.DataFrame, right: pd.DataFrame, columns) -> float:
    distances = []
    for column in columns:
        if _is_categorical(column, left[column], right[column]):
            distances.append(_categorical_distance(left[column], right[column], column))
        else:
            distances.append(_numeric_distance(left[column], right[column]))
    return float(np.mean(distances)) if distances else 0.0


def resample_cat(
    initial: Any,
    regenerated: Mapping[str, Any],
    groups: Sequence[Mapping[str, Any]],
    *,
    columns: Sequence[str],
    seed: int = 0,
) -> Dict[str, Any]:
    initial_cat = _frame(initial)
    if "sample" not in initial_cat.columns:
        raise ValueError("initial CAT needs a sample column")
    known = set(initial_cat["sample"].astype(str))
    pools = {str(name): _frame(table) for name, table in regenerated.items()}
    selected = _validate_columns(initial_cat, pools, columns)
    rng = np.random.default_rng(int(seed))
    parts = []
    reports = []
    names = [str(group.get("name", "")) for group in groups]
    if not names or any(not name for name in names):
        raise ValueError("resampling groups need non-empty names")
    if len(names) != len(set(names)):
        raise ValueError("resampling group names must be unique")
    owner: Dict[str, str] = {}
    for group in groups:
        name = str(group["name"])
        for sample in group.get("samples") or []:
            key = str(sample)
            if key in owner:
                raise ValueError(f"sample {key} is in groups {owner[key]} and {name}")
            owner[key] = name
    for group in groups:
        name = str(group["name"])
        samples = [str(sample) for sample in group["samples"]]
        if not samples:
            raise ValueError(f"group {name} needs at least one initial sample")
        unknown = [sample for sample in samples if sample not in known]
        if unknown:
            raise ValueError(
                f"group {name} sample ids are not in the initial CAT: {', '.join(unknown)}"
            )
        dataset = str(group["dataset"])
        if dataset not in pools:
            raise ValueError(f"group {name} dataset {dataset!r} is not in the regenerated CAT")
        pool = pools[dataset]
        collision = set(PROVENANCE_COLUMNS) & set(pool.columns)
        if collision:
            raise ValueError(
                f"dataset {dataset} already has reserved resampling columns: "
                + ", ".join(sorted(collision))
            )
        member_rows = initial_cat[initial_cat["sample"].astype(str).isin(samples)]
        if "size" not in group:
            raise ValueError(f"group {name} needs an explicit regenerated row size")
        n_rows = int(group["size"])
        if n_rows < 0:
            raise ValueError(f"group {name} size must be non-negative")
        if n_rows and pool.empty:
            raise ValueError(f"dataset {dataset} has no rows to draw")
        descriptors = [_descriptor(column, member_rows[column], pool[column]) for column in selected]
        joint = _joint_descriptor(descriptors)
        weights = _rake_weights([*descriptors, joint], len(pool)) if len(pool) else np.asarray([])
        best = None
        best_mismatch = math.inf
        for _ in range(8 if n_rows else 1):
            candidate = (
                rng.choice(len(pool), size=n_rows, replace=True, p=weights)
                if n_rows
                else np.asarray([], dtype=int)
            )
            mismatch = (
                distribution_mismatch(member_rows, pool.iloc[candidate], selected)
                if n_rows
                else 0.0
            )
            if mismatch < best_mismatch:
                best, best_mismatch = candidate, mismatch
        chosen = best if best is not None else np.asarray([], dtype=int)
        drawn = pool.iloc[list(chosen)].copy().reset_index(drop=True)
        drawn.insert(0, "draw_id", [f"{name}:{i}" for i in range(n_rows)])
        drawn.insert(1, "source_row", [f"{dataset}:{int(i)}" for i in chosen])
        drawn.insert(2, "resampling_group", name)
        parts.append(drawn)
        gaps = [
            {
                "column": descriptor["column"],
                "missing_count": len(descriptor["missing"]),
                "examples": [str(value) for value in descriptor["missing"][:20]],
            }
            for descriptor in [*descriptors, joint]
            if descriptor["missing"]
        ]
        reports.append(
            {
                "group": name,
                "dataset": dataset,
                "n": n_rows,
                "samples": samples,
                "support_gaps": gaps,
                "residual_mismatch": best_mismatch,
            }
        )
    table = pd.concat(parts, ignore_index=True) if parts else initial_cat.iloc[0:0].copy()
    return {
        "schema_version": SCHEMA_VERSION,
        "solver": SOLVER,
        "seed": int(seed),
        "columns": selected,
        "table": table,
        "groups": reports,
        "native_profiler_resampled": False,
    }


def group_distributions(initial, groups, columns) -> Dict[str, Dict[str, Dict[str, float]]]:
    """Frequency of each matching column on one group's initial rows."""
    frame = initial if isinstance(initial, pd.DataFrame) else _initial_frame(initial)
    if "sample" not in frame.columns:
        raise ValueError("initial CAT needs a sample column")
    selected = [str(column) for column in columns]
    forbidden = [column for column in selected if column.strip().lower() in _FORBIDDEN_MATCH]
    if forbidden:
        raise ValueError("resampling matching cannot use " + ", ".join(forbidden))
    missing = [column for column in selected if column not in frame.columns]
    if missing:
        raise ValueError(
            "distribution columns are not in the initial CAT: " + ", ".join(missing)
        )
    distributions: Dict[str, Dict[str, Dict[str, float]]] = {}
    for group in groups:
        name = str(group["name"])
        samples = {str(sample) for sample in group["samples"]}
        rows = frame[frame["sample"].astype(str).isin(samples)]
        column_dist: Dict[str, Dict[str, float]] = {}
        for column in selected:
            counts = rows[column].astype(str).value_counts(dropna=False)
            total = float(counts.sum()) or 1.0
            column_dist[column] = {
                str(value): float(counts[value]) / total for value in counts.index
            }
        distributions[name] = column_dist
    return distributions


def _initial_frame(initial: Mapping[str, pd.DataFrame]) -> pd.DataFrame:
    frames = []
    for sample, table in initial.items():
        work = _frame(table)
        work["sample"] = str(sample)
        frames.append(work)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def save_resampling(result: Mapping[str, Any], output_dir) -> None:
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    table = result["table"]
    table.to_csv(destination / "resampled_training_table.csv", index=False)
    report = {key: value for key, value in result.items() if key != "table"}
    (destination / "resampling.json").write_text(
        json.dumps(report, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    from samovar.abundance import annotation_to_abundance_tables

    abundance_dir = destination / "resampled_training_abundance"
    abundance_dir.mkdir(exist_ok=True)
    counted = table[
        [
            column
            for column in table.columns
            if not str(column).lower().startswith("taxid_samovar")
        ]
    ]
    for name, abundance in annotation_to_abundance_tables(counted).items():
        abundance.to_csv(abundance_dir / f"{name}.csv", index=False)


def _run_resampler(tool, initial, regenerated, groups, config) -> dict:
    if tool in {"builtin", "native", "default"}:
        result = resample_cat(
            initial,
            regenerated,
            groups,
            columns=config.get("columns") or [],
            seed=int(config.get("seed", 0)),
        )
    else:
        from samovar.main_config import iter_tools, parse_tool_entry, tool_path
        from samovar.paths import load_config
        from samovar.tool_contracts import load_python_module

        tools = iter_tools(load_config())
        if tool not in tools:
            raise ValueError(
                f"unknown resampling tool {tool!r}; import it with "
                "`samovar tools import --type resampling`"
            )
        spec = parse_tool_entry(tools[tool], tool)
        if len(spec) < 4 or spec[3] != "resampling":
            raise ValueError(f"{tool!r} is not a resampling tool")
        module = load_python_module(Path(tool_path(spec, tool)), f"samovar_resampling_{tool}")
        function = getattr(module, "resample", None)
        if not callable(function):
            raise ValueError(f"resampling tool {tool!r} does not define resample()")
        result = function(initial, regenerated, groups, dict(config))
    if not isinstance(result, dict) or not isinstance(result.get("table"), pd.DataFrame):
        raise ValueError(f"resampling tool {tool!r} must return a dict with DataFrame table")
    missing = {
        "schema_version",
        "solver",
        "seed",
        "columns",
        "groups",
        "native_profiler_resampled",
    } - set(result)
    if missing:
        raise ValueError(f"resampling tool {tool!r} omitted: {', '.join(sorted(missing))}")
    table = result["table"]
    if not set(PROVENANCE_COLUMNS) <= set(table.columns):
        raise ValueError(f"resampling tool {tool!r} omitted draw provenance columns")
    if not table["draw_id"].is_unique:
        raise ValueError(f"resampling tool {tool!r} returned duplicate draw_id values")
    return result


def apply_resampling(
    regenerated: pd.DataFrame,
    initial: Mapping[str, pd.DataFrame],
    config: Mapping[str, Any],
    output_dir,
) -> tuple[pd.DataFrame, dict | None]:
    block = config.get("resampling")
    if not isinstance(block, Mapping) or not block:
        return regenerated, None
    tool = str(block.get("tool") or "builtin")
    tool_config = dict(block)
    tool_config["seed"] = int(block.get("seed", config.get("seed", 0)))
    result = _run_resampler(
        tool,
        _initial_frame(initial),
        {"regenerated": regenerated.copy()},
        block.get("groups") or [],
        tool_config,
    )
    save_resampling(result, output_dir)
    selected = result["table"].drop(columns=["draw_id", "source_row"], errors="ignore")
    return selected, {key: value for key, value in result.items() if key != "table"}
