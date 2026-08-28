"""ML reprofilers: builtins (RF / AdaBoost / ensemble) plus imported tools.

Train on regenerated annotation tables (known ``true`` taxid) and optional
community ground-truth abundance tables, then profile the initial annotation
tables. Custom tools: ``samovar tools import --type ml``.
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union

import pandas as pd
import yaml

from samovar.main_config import (
    flags_target_matches,
    iter_tools,
    merge_flag_strings,
    parse_tool_entry,
    tool_flags,
    tool_path,
)
from samovar.paths import load_config
from samovar.reprofiling import (
    annotator_taxid_columns,
    load_read_features,
    merge_read_features,
    passthrough_reprofile_tables,
    plot_roc_curves,
    predict_taxid,
    preprocess_data,
    save_model,
    train_models,
)
from samovar.table_regenerators import extra_flags_argv

REPROFILER_GROUP = "reprofiler"
PathLike = Union[str, Path]

_BUILTIN_ALIASES = {
    "ensemble": "ensemble",
    "samovar": "ensemble",
    "default": "ensemble",
    "both": "ensemble",
    "random_forest": "random_forest",
    "rf": "random_forest",
    "randomforest": "random_forest",
    "adaboost": "adaboost",
    "ada": "adaboost",
}

REPROFILER_FLAG_GROUPS = (
    "reprofiler",
    "reprofiling",
    "reprofile",
    "ml",
)


class MissingReprofilerError(ValueError):
    """``tools.<name>`` is missing or is not a reprofiler."""


@dataclass
class ReprofileResult:
    """Tables to write plus the trained model (and optional extras)."""

    tables: Dict[str, pd.DataFrame]
    model: Any = None
    models: Optional[Dict[str, Any]] = None
    metrics: Optional[Dict[str, float]] = None
    feature_cols: Optional[List[str]] = None
    extras: Dict[str, Any] = field(default_factory=dict)


def resolve_reprofiler(name: Optional[str]) -> Tuple[str, str]:
    key = str(name or "ensemble").strip() or "ensemble"
    low = key.lower().replace("-", "_")
    if low in _BUILTIN_ALIASES:
        return "builtin", _BUILTIN_ALIASES[low]
    return "custom", key


def require_known_reprofiler(name: Optional[str]) -> str:
    kind, canon = resolve_reprofiler(name)
    if kind == "custom":
        lookup_reprofiler(canon)
        return canon
    return canon


def flags_apply_to_reprofiler(target: str, name: Optional[str] = None) -> bool:
    try:
        kind, canon = resolve_reprofiler(name)
    except MissingReprofilerError:
        kind, canon = "custom", name
    names = [canon, name]
    if kind == "builtin":
        if canon == "random_forest":
            names.extend(["rf", "randomforest"])
        if canon == "adaboost":
            names.extend(["ada"])
        if canon == "ensemble":
            names.extend(["samovar", "default", "both"])
    return flags_target_matches(target, *names, groups=REPROFILER_FLAG_GROUPS)


def is_reprofiler_flag_target(target: str) -> bool:
    extra = list(_BUILTIN_ALIASES.keys())
    extra.extend(iter_custom_reprofiler_names())
    return flags_target_matches(target, *extra, groups=REPROFILER_FLAG_GROUPS)


def iter_custom_reprofiler_names(tools: Optional[Dict[str, List[str]]] = None) -> List[str]:
    mapping = tools if tools is not None else iter_tools(load_config())
    names: List[str] = []
    for name, spec in mapping.items():
        parsed = parse_tool_entry(spec, name)
        if str(parsed[3] or "").strip() != REPROFILER_GROUP:
            continue
        if name.lower() in _BUILTIN_ALIASES:
            continue
        if tool_path(parsed, name):
            names.append(name)
    return names


def lookup_reprofiler(name: str) -> list:
    key = str(name or "").strip()
    if not key:
        raise MissingReprofilerError(
            "Empty reprofiler name. Import with "
            "`samovar tools import -n NAME --type ml` or use ensemble|random_forest|adaboost."
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
        raise MissingReprofilerError(
            f"reprofiler {key!r} is not in the install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type ml`."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != REPROFILER_GROUP:
        raise MissingReprofilerError(
            f"tools.{matched} has group {group!r}, expected {REPROFILER_GROUP!r}. "
            "Re-import with --type ml (or --type reprofiler)."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingReprofilerError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return parsed


def attach_reprofiler_flags(name: str, config: Dict[str, Any]) -> Dict[str, Any]:
    cfg = dict(config or {})
    kind, canon = resolve_reprofiler(name)
    imported = ""
    if kind == "custom":
        spec = lookup_reprofiler(canon)
        imported = tool_flags(spec, canon)
    named_map = cfg.get("reprofiler_tool_flags") or {}
    named = named_map.get(canon) or named_map.get(name)
    if named is None:
        low = str(canon).lower()
        for key, value in named_map.items():
            if str(key).lower() == low:
                named = value
                break
    cfg["extra_flags"] = merge_flag_strings(
        imported,
        cfg.get("reprofiler_flags"),
        named,
        cfg.get("extra_flags"),
    )
    cfg["extra_argv"] = extra_flags_argv(cfg["extra_flags"])
    return cfg


def load_csv_tables(
    directory: Optional[PathLike],
    *,
    skip_prefixes: Sequence[str] = ("combined_annotation_table",),
) -> Dict[str, pd.DataFrame]:
    if not directory:
        return {}
    root = Path(directory)
    if not root.is_dir():
        return {}
    skip = tuple(skip_prefixes)
    tables: Dict[str, pd.DataFrame] = {}
    for path in sorted(root.glob("*.csv")):
        if any(path.name.startswith(prefix) for prefix in skip):
            continue
        tables[path.stem] = pd.read_csv(path)
    return tables


def load_regenerated_table(path: Optional[PathLike]) -> pd.DataFrame:
    if not path:
        return pd.DataFrame()
    target = Path(path)
    if target.is_file():
        return pd.read_csv(target)
    if target.is_dir():
        combined = target / "combined_annotation_table.csv"
        if combined.is_file():
            return pd.read_csv(combined)
        parts = list(load_csv_tables(target).values())
        if parts:
            return pd.concat(parts, ignore_index=True)
    return pd.DataFrame()


def profile_tables_with_model(
    initial: Dict[str, pd.DataFrame],
    model: Any,
    *,
    features_df=None,
    classify_unclassified: bool = False,
) -> Dict[str, pd.DataFrame]:
    out: Dict[str, pd.DataFrame] = {}
    for stem, raw in initial.items():
        df = merge_read_features(raw.copy(), features_df)
        fill_cols = [c for c in df.columns if c != "true"]
        if fill_cols:
            df[fill_cols] = df[fill_cols].fillna(0)
        result = predict_taxid(df, model_path=model)
        if not classify_unclassified:
            taxid_cols = annotator_taxid_columns(df)
            if taxid_cols:
                numeric = df[taxid_cols].apply(pd.to_numeric, errors="coerce").fillna(0)
                mask = (numeric == 0).all(axis=1)
                result.loc[mask, "taxid_SAMOVAR"] = 0
                result.loc[mask, "taxid_SAMOVAR_confidence"] = 0
        out[_reprofiled_stem(stem)] = result
    return out


def _reprofiled_stem(stem: str) -> str:
    root = Path(stem).name.split(".")[0]
    if root.endswith("_reprofiled"):
        return root
    return f"{root}_reprofiled"


def write_reprofile_result(result: ReprofileResult, output_dir: PathLike) -> Path:
    dest = Path(output_dir)
    dest.mkdir(parents=True, exist_ok=True)
    if result.model is not None:
        save_model(result.model, str(dest / "trained_model.joblib"))
    for stem, table in (result.tables or {}).items():
        name = f"{_reprofiled_stem(stem)}.csv"
        table.to_csv(dest / name, index=False)
    return dest


def _coerce_result(value: Any, initial: Dict[str, pd.DataFrame]) -> ReprofileResult:
    if isinstance(value, ReprofileResult):
        return value
    if isinstance(value, dict):
        tables = value.get("tables") or value.get("profiles") or {}
        return ReprofileResult(
            tables=tables,
            model=value.get("model"),
            models=value.get("models"),
            metrics=value.get("metrics"),
            feature_cols=value.get("feature_cols"),
            extras={k: v for k, v in value.items() if k not in {
                "tables", "profiles", "model", "models", "metrics", "feature_cols"
            }},
        )
    raise TypeError(
        f"reprofile() must return ReprofileResult or dict, got {type(value).__name__}"
    )


def run_builtin_reprofiler(
    name: str,
    regenerated: pd.DataFrame,
    ground_truth: Dict[str, pd.DataFrame],
    initial: Dict[str, pd.DataFrame],
    config: Dict[str, Any],
) -> ReprofileResult:
    _ = ground_truth
    cfg = dict(config or {})
    features_df = cfg.get("features_df")
    validation = merge_read_features(regenerated.copy(), features_df)
    if "seq" in validation.columns:
        validation = validation.drop(columns=["seq"])
    if "true" in validation.columns:
        validation = validation.dropna(subset=["true"])
    fill_cols = [c for c in validation.columns if c != "true"]
    if fill_cols:
        validation[fill_cols] = validation[fill_cols].fillna(0)
    output_dir = Path(cfg.get("output_dir") or ".")
    output_dir.mkdir(parents=True, exist_ok=True)
    if validation.empty or "true" not in validation.columns:
        passthrough_reprofile_tables(
            str(cfg.get("initial_dir") or output_dir),
            str(output_dir),
            features_df=features_df,
        )
        return ReprofileResult(tables=load_csv_tables(output_dir, skip_prefixes=()))
    seed = int(cfg.get("seed") or 42)
    try:
        best, models, metrics, feature_cols = train_models(
            validation, random_state=seed, methods=name
        )
    except ValueError:
        passthrough_reprofile_tables(
            str(cfg.get("initial_dir") or output_dir),
            str(output_dir),
            features_df=features_df,
        )
        return ReprofileResult(tables=load_csv_tables(output_dir, skip_prefixes=()))
    tables = profile_tables_with_model(
        initial,
        best,
        features_df=features_df,
        classify_unclassified=bool(cfg.get("classify_unclassified")),
    )
    processed = preprocess_data(validation.copy())
    try:
        X = processed.reindex(columns=feature_cols, fill_value=0)
        plot_roc_curves(models, X, processed["true"], output_dir=str(output_dir))
    except Exception:
        pass
    return ReprofileResult(
        tables=tables,
        model=best,
        models=models,
        metrics=metrics,
        feature_cols=feature_cols,
    )


def run_reprofiler(
    name: Optional[str],
    *,
    regenerated: Optional[pd.DataFrame] = None,
    regenerated_path: Optional[PathLike] = None,
    ground_truth: Optional[Dict[str, pd.DataFrame]] = None,
    ground_truth_dir: Optional[PathLike] = None,
    initial: Optional[Dict[str, pd.DataFrame]] = None,
    initial_dir: Optional[PathLike] = None,
    output_dir: PathLike,
    config: Optional[Dict[str, Any]] = None,
) -> ReprofileResult:
    kind, canon = resolve_reprofiler(name)
    cfg = attach_reprofiler_flags(canon if kind == "custom" else canon, dict(config or {}))
    cfg["output_dir"] = str(output_dir)
    cfg["seed"] = int(cfg.get("seed") or 42)
    if initial_dir:
        cfg["initial_dir"] = str(initial_dir)
    if regenerated_path:
        cfg["regenerated_path"] = str(regenerated_path)
    if ground_truth_dir:
        cfg["ground_truth_dir"] = str(ground_truth_dir)

    regen = regenerated if regenerated is not None else load_regenerated_table(
        regenerated_path
    )
    truth = ground_truth if ground_truth is not None else load_csv_tables(
        ground_truth_dir, skip_prefixes=()
    )
    init = initial if initial is not None else load_csv_tables(initial_dir)
    cfg["features_df"] = cfg.get("features_df")
    if cfg.get("features") and cfg["features_df"] is None:
        cfg["features_df"] = load_read_features(cfg.get("features"))

    if kind == "builtin":
        result = run_builtin_reprofiler(canon, regen, truth, init, cfg)
    else:
        result = _run_custom(canon, regen, truth, init, cfg)
    if not result.tables:
        src = initial_dir or cfg.get("initial_dir")
        if src:
            passthrough_reprofile_tables(
                str(src),
                str(output_dir),
                features_df=cfg.get("features_df"),
            )
        result.tables = load_csv_tables(output_dir, skip_prefixes=())
        return result
    write_reprofile_result(result, output_dir)
    return result


def _run_custom(
    name: str,
    regenerated: pd.DataFrame,
    ground_truth: Dict[str, pd.DataFrame],
    initial: Dict[str, pd.DataFrame],
    config: Dict[str, Any],
) -> ReprofileResult:
    spec = lookup_reprofiler(name)
    path = Path(tool_path(spec, name)).expanduser()
    py_fn = _try_python_callable(path, name)
    if py_fn is not None:
        raw = py_fn(regenerated, ground_truth, initial, config)
        return _coerce_result(raw, initial)
    return _run_cli_reprofiler(path, regenerated, ground_truth, initial, config)


def _try_python_callable(path: Path, name: str) -> Optional[Callable]:
    if path.suffix.lower() != ".py" or not path.is_file():
        return None
    try:
        spec = importlib.util.spec_from_file_location(
            f"samovar_custom_reprofiler_{name}", path
        )
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except Exception:
        return None
    fn = getattr(module, "reprofile", None)
    if callable(fn):
        return fn
    cls = getattr(module, "Reprofiler", None) or getattr(module, "ReprofilerModel", None)
    if cls is None:
        return None
    try:
        inst = cls()
    except Exception:
        return None
    run = getattr(inst, "reprofile", None) or getattr(inst, "run", None)
    return run if callable(run) else None


def _run_cli_reprofiler(
    path: Path,
    regenerated: pd.DataFrame,
    ground_truth: Dict[str, pd.DataFrame],
    initial: Dict[str, pd.DataFrame],
    config: Dict[str, Any],
) -> ReprofileResult:
    output_dir = Path(config.get("output_dir") or ".")
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp = output_dir / ".reprofiler_inputs"
    tmp.mkdir(parents=True, exist_ok=True)
    regen_path = tmp / "regenerated.csv"
    gt_dir = tmp / "ground_truth"
    init_dir = tmp / "initial"
    gt_dir.mkdir(exist_ok=True)
    init_dir.mkdir(exist_ok=True)
    regenerated.to_csv(regen_path, index=False)
    for stem, table in ground_truth.items():
        table.to_csv(gt_dir / f"{stem}.csv", index=False)
    for stem, table in initial.items():
        table.to_csv(init_dir / f"{stem}.csv", index=False)
    cmd = [
        str(path),
        "--regenerated",
        str(regen_path),
        "--ground-truth",
        str(gt_dir),
        "--initial",
        str(init_dir),
        "-o",
        str(output_dir),
    ]
    extra = list(config.get("extra_argv") or extra_flags_argv(config.get("extra_flags")))
    cmd.extend(extra)
    if path.suffix.lower() == ".py":
        cmd = [sys.executable, *cmd]
    subprocess.run(cmd, check=True)
    tables = load_csv_tables(output_dir, skip_prefixes=())
    model_path = output_dir / "trained_model.joblib"
    model = str(model_path) if model_path.is_file() else None
    return ReprofileResult(tables=tables, model=model)


def load_reprofiler_config(path: Optional[PathLike]) -> Dict[str, Any]:
    if not path:
        return {}
    file = Path(path)
    if not file.is_file():
        return {}
    data = yaml.safe_load(file.read_text(encoding="utf-8")) or {}
    return data if isinstance(data, dict) else {}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="SamovaR ML reprofiling")
    parser.add_argument(
        "--reprofiling_dir",
        "-r",
        dest="initial_dir",
        help="Initial annotation tables to profile",
    )
    parser.add_argument(
        "--initial",
        dest="initial_dir_alt",
        default=None,
        help="Alias of --reprofiling_dir",
    )
    parser.add_argument(
        "--validation_file",
        "-v",
        dest="regenerated",
        help="Regenerated annotation table (known true taxid)",
    )
    parser.add_argument(
        "--regenerated",
        dest="regenerated_alt",
        default=None,
        help="Alias of --validation_file",
    )
    parser.add_argument(
        "--ground-truth",
        "--ground_truth",
        dest="ground_truth",
        default=None,
        help="Directory of abundance CSVs used to generate regenerated metagenomes",
    )
    parser.add_argument("--output_dir", "-o", required=True)
    parser.add_argument(
        "--reprofiler",
        "--ml",
        dest="reprofiler",
        default=None,
        help="ensemble|random_forest|adaboost or an imported --type ml name",
    )
    parser.add_argument("--config", default="", help="config_reprofiling.yaml from prepare")
    parser.add_argument(
        "--classify-unclassified",
        action="store_true",
        help="Classify reads where every annotator taxid is 0",
    )
    parser.add_argument("--features", "-f", default=None)
    parser.add_argument("--seed", type=int, default=None)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    cfg = load_reprofiler_config(args.config)
    initial_dir = args.initial_dir or args.initial_dir_alt or cfg.get("initial_dir")
    regenerated = args.regenerated or args.regenerated_alt or cfg.get("regenerated_path")
    if not initial_dir:
        print("Error: --reprofiling_dir / --initial is required", file=sys.stderr)
        return 1
    if not regenerated:
        print("Error: --validation_file / --regenerated is required", file=sys.stderr)
        return 1
    name = args.reprofiler or cfg.get("reprofiler") or "ensemble"
    if args.seed is not None:
        cfg["seed"] = args.seed
    if args.classify_unclassified:
        cfg["classify_unclassified"] = True
    if args.features:
        cfg["features"] = args.features
    cfg.setdefault("seed", 42)
    try:
        require_known_reprofiler(name)
        run_reprofiler(
            name,
            regenerated_path=regenerated,
            ground_truth_dir=args.ground_truth or cfg.get("ground_truth_dir"),
            initial_dir=initial_dir,
            output_dir=args.output_dir,
            config=cfg,
        )
    except (MissingReprofilerError, subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
