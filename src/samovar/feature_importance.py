"""Feature-importance scoring at reprofiler training time.

Universal extractor for sklearn-like estimators and later analogues (HMM,
GMM, any object with ``predict`` / ``coef_`` / ``emissionprob_``).

Imported tools: ``samovar tools import --type feature-importance`` with
``score_feature_importance(model, annotation, initial_abundance,
regenerated_abundance, config)``.
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from samovar.main_config import (
    flags_target_matches,
    iter_tools,
    parse_tool_entry,
    tool_flags,
    tool_path,
)
from samovar.paths import load_config

PathLike = Union[str, Path]
FEATURE_IMPORTANCE_GROUP = "feature_importance"
NONE_TOKENS = frozenset({"", "none", "off", "false", "0", "skip", "no"})
BUILTIN_NAMES = frozenset(
    {
        "builtin",
        "native",
        "sklearn",
        "permutation",
        "hmm",
        "default",
        "auto",
    }
)
FEATURE_IMPORTANCE_FLAG_GROUPS = (
    "feature_importance",
    "feature-importance",
    "featureimportance",
    "fi",
    "importance",
)


class MissingFeatureImportanceError(ValueError):
    """``tools.<name>`` is missing or is not feature-importance scoring."""


def flags_apply_to_feature_importance(target: str, name: Optional[str] = None) -> bool:
    names = [name] if name else []
    return flags_target_matches(
        target,
        *names,
        groups=FEATURE_IMPORTANCE_FLAG_GROUPS,
    )


def is_feature_importance_flag_target(target: str) -> bool:
    extra = list(iter_custom_feature_importance_names())
    extra.extend(BUILTIN_NAMES)
    extra.append("builtin")
    return flags_target_matches(target, *extra, groups=FEATURE_IMPORTANCE_FLAG_GROUPS)


def iter_custom_feature_importance_names(tools: Optional[Dict[str, list]] = None) -> List[str]:
    names: List[str] = []
    mapping = tools if tools is not None else iter_tools(load_config())
    for key, spec in mapping.items():
        parsed = parse_tool_entry(spec, key)
        if len(parsed) > 3 and parsed[3] == FEATURE_IMPORTANCE_GROUP:
            names.append(key)
    return names


def canonicalize_feature_importance(name: Optional[str]) -> str:
    key = str(name or "builtin").strip()
    low = key.lower().replace("-", "_")
    if low in NONE_TOKENS:
        return "none"
    if low in BUILTIN_NAMES:
        return "builtin"
    try:
        matched, _spec = lookup_feature_importance(key)
        return matched
    except MissingFeatureImportanceError:
        pass
    raise MissingFeatureImportanceError(
        f"Unknown feature-importance scorer {name!r}. Built-in: builtin. "
        "Or import one with `samovar tools import --type feature-importance`."
    )


def lookup_feature_importance(name: str) -> Tuple[str, list]:
    from samovar.baselines import BASELINE_TOOLS

    key = str(name or "").strip()
    if not key:
        raise MissingFeatureImportanceError("Empty feature-importance scorer name.")
    low = key.lower().replace("-", "_")
    if low in BUILTIN_NAMES:
        path = BASELINE_TOOLS[FEATURE_IMPORTANCE_GROUP]
        return "builtin", ["", "bash", str(path), FEATURE_IMPORTANCE_GROUP]
    tools = iter_tools(load_config())
    spec = tools.get(key)
    matched = key
    if spec is None:
        for stored, row in tools.items():
            if stored.lower().replace("-", "_") == low:
                spec = row
                matched = stored
                break
    if spec is None:
        raise MissingFeatureImportanceError(
            f"feature-importance {key!r} is not in the install config. "
            "Import with `samovar tools import --type feature-importance`."
        )
    parsed = parse_tool_entry(spec, matched)
    group = parsed[3] if len(parsed) > 3 else ""
    if group and group != FEATURE_IMPORTANCE_GROUP:
        raise MissingFeatureImportanceError(
            f"{matched!r} is type {group!r}, not feature_importance. "
            "Re-import with --type feature-importance."
        )
    return matched, parsed


def unwrap_estimator(model: Any) -> Any:
    """Walk wrappers (joblib path, Pipeline, CalibratedClassifier, OvR)."""
    current = model
    if isinstance(current, (str, Path)):
        path = Path(current)
        if path.is_file():
            import joblib

            current = joblib.load(path)
    seen = set()
    for _ in range(8):
        ident = id(current)
        if ident in seen or current is None:
            break
        seen.add(ident)
        if _has_native_importance(current):
            break
        if hasattr(current, "steps") and current.steps:
            current = current.steps[-1][1]
            continue
        if hasattr(current, "best_estimator_") and current.best_estimator_ is not None:
            current = current.best_estimator_
            continue
        cals = getattr(current, "calibrated_classifiers_", None)
        if cals:
            first = cals[0]
            current = getattr(first, "estimator", None) or getattr(first, "base_estimator", first)
            continue
        nested = getattr(current, "estimator_", None) or getattr(current, "base_estimator_", None)
        if nested is not None and nested is not current:
            current = nested
            continue
        break
    return current


def _has_native_importance(est: Any) -> bool:
    return any(
        hasattr(est, attr)
        for attr in (
            "feature_importances_",
            "feature_importance_",
            "importances_",
            "coef_",
            "emissionprob_",
            "emissprob_",
            "means_",
        )
    )


def _as_1d(values: Any, n_features: Optional[int] = None) -> Optional[np.ndarray]:
    if values is None:
        return None
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return None
    if arr.ndim > 1:
        arr = np.abs(arr).mean(axis=tuple(range(arr.ndim - 1)))
    arr = np.ravel(arr)
    if n_features is not None and arr.size != n_features:
        return None
    arr = np.where(np.isfinite(arr), arr, 0.0)
    return arr


def extract_native_importance(
    model: Any,
    feature_names: Optional[Sequence[str]] = None,
) -> Tuple[Optional[np.ndarray], str]:
    """Pull a native importance vector from an arbitrary fitted object."""
    est = unwrap_estimator(model)
    names = list(feature_names or getattr(est, "feature_names_in_", None) or [])
    n_feat = len(names) if names else None

    getter = getattr(est, "get_feature_importance", None) or getattr(
        est, "get_feature_importances", None
    )
    if callable(getter):
        try:
            raw = getter()
            vec = _as_1d(raw, n_feat)
            if vec is not None:
                return vec, "get_feature_importance"
        except Exception:
            pass

    for attr in (
        "feature_importances_",
        "feature_importance_",
        "importances_",
        "variable_importances_",
    ):
        if hasattr(est, attr):
            vec = _as_1d(getattr(est, attr), n_feat)
            if vec is not None:
                return vec, attr

    if hasattr(est, "coef_"):
        vec = _as_1d(getattr(est, "coef_"), n_feat)
        if vec is not None:
            return vec, "coef_"

    for attr in ("emissionprob_", "emissprob_", "means_"):
        if hasattr(est, attr):
            raw = np.asarray(getattr(est, attr), dtype=float)
            if raw.ndim >= 2:
                vec = _as_1d(np.var(raw, axis=0), n_feat)
            else:
                vec = _as_1d(raw, n_feat)
            if vec is not None:
                return vec, attr

    trans = getattr(est, "transmat_", None)
    if trans is not None and n_feat:
        # HMM with no per-feature emissions: spread transition mass equally.
        mass = float(np.abs(np.asarray(trans, dtype=float)).sum())
        return np.full(n_feat, mass / max(n_feat, 1), dtype=float), "transmat_"

    return None, ""


def permutation_importance_scores(
    model: Any,
    X: pd.DataFrame,
    y: pd.Series,
    *,
    n_repeats: int = 5,
    random_state: int = 0,
    sample_weight: Optional[np.ndarray] = None,
) -> Optional[np.ndarray]:
    if X is None or X.empty or y is None or len(y) != len(X):
        return None
    est = unwrap_estimator(model)
    predict = getattr(est, "predict", None)
    if not callable(predict):
        return None
    try:
        from sklearn.inspection import permutation_importance
        from sklearn.metrics import accuracy_score

        def _acc(estimator, features, labels):
            pred = estimator.predict(features)
            if sample_weight is None:
                return float(accuracy_score(labels, pred))
            return float(accuracy_score(labels, pred, sample_weight=sample_weight))

        result = permutation_importance(
            est,
            X,
            y,
            n_repeats=max(1, int(n_repeats)),
            random_state=int(random_state),
            scoring=_acc,
        )
        return _as_1d(result.importances_mean, X.shape[1])
    except Exception:
        pass
    rng = np.random.default_rng(int(random_state))
    try:
        pred = np.asarray(predict(X))
    except Exception:
        return None
    labels = np.asarray(y)
    weights = sample_weight if sample_weight is not None else np.ones(len(labels), dtype=float)
    weights = np.asarray(weights, dtype=float)
    if weights.shape[0] != len(labels):
        weights = np.ones(len(labels), dtype=float)

    def _weighted_acc(pred_labels) -> float:
        return float(np.average((np.asarray(pred_labels) == labels).astype(float), weights=weights))

    baseline = _weighted_acc(pred)
    drops = np.zeros(X.shape[1], dtype=float)
    work = X.copy()
    repeats = max(1, int(n_repeats))
    for j, col in enumerate(X.columns):
        accs = []
        original = work[col].to_numpy().copy()
        for _ in range(repeats):
            work[col] = rng.permutation(original)
            try:
                accs.append(_weighted_acc(predict(work)))
            except Exception:
                accs.append(baseline)
        work[col] = original
        drops[j] = max(0.0, baseline - float(np.mean(accs)))
    return drops


def _as_frame(value: Any) -> pd.DataFrame:
    if value is None:
        return pd.DataFrame()
    if isinstance(value, pd.DataFrame):
        return value
    to_long = getattr(value, "to_long_table", None) or getattr(value, "to_dataframe", None)
    if callable(to_long):
        frame = to_long()
        if isinstance(frame, pd.DataFrame):
            return frame
    table = getattr(value, "table", None)
    if isinstance(table, pd.DataFrame):
        return table
    if isinstance(value, (str, Path)) and Path(value).is_file():
        path = Path(value)
        sep = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
        return pd.read_csv(path, sep=sep)
    return pd.DataFrame()


def _as_tables(value: Any) -> Dict[str, pd.DataFrame]:
    if value is None:
        return {}
    if isinstance(value, dict):
        out = {}
        for key, table in value.items():
            frame = _as_frame(table)
            if not frame.empty:
                out[str(key)] = frame
        return out
    path = Path(value) if isinstance(value, (str, Path)) else None
    if path is None:
        frame = _as_frame(value)
        return {"table": frame} if not frame.empty else {}
    if path.is_file():
        return {path.stem: _as_frame(path)}
    if path.is_dir():
        tables: Dict[str, pd.DataFrame] = {}
        for child in sorted(path.glob("*.csv")):
            tables[child.stem] = _as_frame(child)
        return tables
    return {}


def _combine_abundance(tables: Dict[str, pd.DataFrame]) -> pd.Series:
    totals: Dict[str, float] = {}
    for table in tables.values():
        if table is None or table.empty:
            continue
        tax_col = "taxid" if "taxid" in table.columns else (
            "taxID" if "taxID" in table.columns else None
        )
        if tax_col is None:
            continue
        sample_cols = [c for c in table.columns if str(c).startswith("N_")]
        if not sample_cols:
            sample_cols = [
                c for c in table.columns if c != tax_col and pd.api.types.is_numeric_dtype(table[c])
            ]
        if not sample_cols:
            continue
        grouped = table.groupby(table[tax_col].astype(str), sort=False)[sample_cols].sum()
        for taxid, row in grouped.iterrows():
            totals[str(taxid)] = totals.get(str(taxid), 0.0) + float(np.nansum(row.to_numpy()))
    series = pd.Series(totals, dtype=float)
    if series.empty or float(series.sum()) <= 0:
        return series
    return series / float(series.sum())


def _predicted_mix(annotation: pd.DataFrame, column: str) -> pd.Series:
    if column not in annotation.columns:
        return pd.Series(dtype=float)
    counts = annotation[column].fillna(0).astype(str).value_counts()
    counts = counts[~counts.index.isin({"0", "nan", "None", ""})]
    if counts.empty or float(counts.sum()) <= 0:
        return pd.Series(dtype=float)
    return counts.astype(float) / float(counts.sum())


def _l1(a: pd.Series, b: pd.Series) -> float:
    idx = sorted(set(a.index) | set(b.index))
    left = a.reindex(idx).fillna(0.0)
    right = b.reindex(idx).fillna(0.0)
    return float(np.abs(left - right).sum()) / 2.0


def profile_agreement_scores(
    annotation: pd.DataFrame,
    feature_names: Sequence[str],
    regenerated: Dict[str, pd.DataFrame],
    initial: Dict[str, pd.DataFrame],
) -> Tuple[Dict[str, float], Dict[str, float]]:
    regen_mix = _combine_abundance(regenerated)
    init_mix = _combine_abundance(initial)
    regen_scores: Dict[str, float] = {}
    init_scores: Dict[str, float] = {}
    for col in feature_names:
        mix = _predicted_mix(annotation, col)
        if mix.empty:
            regen_scores[col] = 0.0
            init_scores[col] = 0.0
            continue
        regen_scores[col] = max(0.0, 1.0 - _l1(mix, regen_mix)) if not regen_mix.empty else 0.0
        init_scores[col] = max(0.0, 1.0 - _l1(mix, init_mix)) if not init_mix.empty else 0.0
    return regen_scores, init_scores


def _sample_weights(annotation: pd.DataFrame, regenerated: Dict[str, pd.DataFrame]) -> Optional[np.ndarray]:
    mix = _combine_abundance(regenerated)
    if mix.empty or "true" not in annotation.columns:
        return None
    mapped = annotation["true"].astype(str).map(mix).fillna(mix.min() if not mix.empty else 0.0)
    weights = mapped.to_numpy(dtype=float)
    if not np.isfinite(weights).any() or float(weights.sum()) <= 0:
        return None
    return weights


def _normalize(vec: np.ndarray) -> np.ndarray:
    total = float(np.abs(vec).sum())
    if total <= 0:
        return np.zeros_like(vec, dtype=float)
    return np.abs(vec) / total


def importance_table(
    *,
    model: Any,
    annotation: Any,
    initial_abundance: Any,
    regenerated_abundance: Any,
    config: Optional[Dict[str, Any]] = None,
    feature_names: Optional[Sequence[str]] = None,
) -> pd.DataFrame:
    cfg = dict(config or {})
    frame = _as_frame(annotation)
    initial = _as_tables(initial_abundance)
    regenerated = _as_tables(regenerated_abundance)
    processed = frame
    y = None
    X = pd.DataFrame()
    try:
        from samovar.reprofiling import preprocess_data

        if not frame.empty and "true" in frame.columns:
            processed = preprocess_data(frame.copy())
            y = processed["true"] if "true" in processed.columns else None
            cols = [
                c
                for c in (feature_names or processed.columns)
                if c in processed.columns and c not in {"seq", "true"}
            ]
            X = processed.reindex(columns=cols, fill_value=0)
    except Exception:
        cols = [c for c in frame.columns if c not in {"seq", "true", "sample"}]
        X = frame.reindex(columns=cols, fill_value=0) if cols else pd.DataFrame()
        y = frame["true"] if "true" in frame.columns else None

    names: List[str] = []
    if feature_names:
        names = [str(c) for c in feature_names]
    elif not X.empty:
        names = [str(c) for c in X.columns]
    else:
        inferred = getattr(unwrap_estimator(model), "feature_names_in_", None)
        if inferred is not None:
            names = [str(c) for c in list(inferred)]
    native, native_src = extract_native_importance(model, names)
    if native is not None and not names:
        names = [f"feature_{i}" for i in range(len(native))]
    if not names:
        return pd.DataFrame(columns=["feature", "importance", "native", "permutation", "source"])

    n_repeats = int(cfg.get("n_repeats") or 5)
    extra = list(cfg.get("extra_argv") or [])
    if "--n-repeats" in extra:
        idx = extra.index("--n-repeats")
        if idx + 1 < len(extra):
            n_repeats = int(extra[idx + 1])
    seed = int(cfg.get("seed") or 0)
    weights = None
    if y is not None and not X.empty:
        aligned = processed.reindex(X.index) if not processed.empty else frame
        weights = _sample_weights(aligned if "true" in aligned.columns else frame, regenerated)
        if weights is not None and len(weights) != len(X):
            weights = None
    perm = None
    force_perm = str(cfg.get("method") or cfg.get("feature_importance") or "").lower() in {
        "permutation"
    } or bool(cfg.get("permutation"))
    native_ok = native is not None and float(np.abs(native).sum()) > 0
    if (force_perm or not native_ok) and not X.empty and y is not None:
        use_x = X.reindex(columns=names, fill_value=0)
        perm = permutation_importance_scores(
            model,
            use_x,
            y,
            n_repeats=n_repeats,
            random_state=seed,
            sample_weight=weights,
        )
    regen_agree, init_agree = profile_agreement_scores(frame, names, regenerated, initial)
    native_map = {}
    if native is not None:
        if len(native) == len(names):
            native_map = dict(zip(names, _normalize(native)))
        elif names:
            native_map = {n: 0.0 for n in names}
    perm_map = {}
    if perm is not None and len(perm) == len(names):
        perm_map = dict(zip(names, _normalize(perm)))

    rows = []
    for name in names:
        native_v = float(native_map.get(name, 0.0))
        perm_v = float(perm_map.get(name, 0.0))
        if native_v > 0:
            importance = native_v
            source = native_src or "native"
        elif perm_v > 0:
            importance = perm_v
            source = "permutation"
        else:
            importance = float(regen_agree.get(name, 0.0))
            source = "profile_agreement" if importance else "none"
        rows.append(
            {
                "feature": name,
                "importance": importance,
                "native": native_v,
                "permutation": perm_v,
                "regen_agreement": float(regen_agree.get(name, 0.0)),
                "initial_agreement": float(init_agree.get(name, 0.0)),
                "source": source,
            }
        )
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    table = table.sort_values("importance", ascending=False, kind="mergesort").reset_index(drop=True)
    return table


def _plot_dir_from_config(config: Dict[str, Any]) -> Path:
    explicit = config.get("plot_dir") or config.get("feature_importance_dir")
    if explicit:
        return Path(explicit)
    output = Path(config.get("output_dir") or ".")
    if output.name == "reprofiled_annotations":
        return output.parent / "feature_importance_plots"
    return output / "feature_importance_plots"


def _write_optional_plots(table: pd.DataFrame, dest: Path, model_name: str) -> List[str]:
    written: List[str] = []
    if table.empty:
        return written
    plot_df = table.head(40).copy()
    html_path = dest / f"feature_importance_{_slug(model_name)}.html"
    png_path = dest / f"feature_importance_{_slug(model_name)}.png"
    svg_path = dest / f"feature_importance_{_slug(model_name)}.svg"
    try:
        import altair as alt

        chart = (
            alt.Chart(plot_df)
            .mark_bar()
            .encode(
                x=alt.X("importance:Q", title="Importance"),
                y=alt.Y("feature:N", sort="-x", title="Feature"),
                tooltip=["feature", "importance", "native", "permutation", "source"],
            )
            .properties(title=f"Feature importance ({model_name})", height=max(120, 16 * len(plot_df)))
        )
        chart.save(str(html_path))
        written.append(str(html_path))
    except Exception:
        pass
    try:
        import cnsplots as cns

        cns.figure(width=420, height=max(180, 18 * len(plot_df)))
        bar = getattr(cns, "barplot", None) or getattr(cns, "bar", None)
        if callable(bar):
            try:
                bar(data=plot_df, x="importance", y="feature")
            except TypeError:
                bar(plot_df["importance"], plot_df["feature"])
        else:
            import matplotlib.pyplot as plt

            plt.barh(plot_df["feature"][::-1], plot_df["importance"][::-1])
        save = getattr(cns, "savefig", None)
        if callable(save):
            save(str(png_path))
            try:
                save(str(svg_path))
                written.append(str(svg_path))
            except Exception:
                pass
        written.append(str(png_path))
    except Exception:
        pass
    return written


def _slug(text: str) -> str:
    keep = []
    for ch in str(text):
        if ch.isalnum() or ch in {"-", "_"}:
            keep.append(ch)
        else:
            keep.append("_")
    return "".join(keep).strip("_") or "model"


def write_feature_importance_outputs(
    table: pd.DataFrame,
    dest: PathLike,
    *,
    model_name: str = "best",
    description: str = "",
) -> Dict[str, Any]:
    from samovar.stage_report import write_bargraph_mqc, write_table_mqc

    folder = Path(dest)
    folder.mkdir(parents=True, exist_ok=True)
    model_slug = _slug(model_name)
    tsv = folder / f"feature_importance_{model_slug}.tsv"
    table.to_csv(tsv, sep="\t", index=False)
    rows = table.to_dict(orient="records")
    write_table_mqc(
        rows,
        folder / f"feature_importance_{model_slug}_mqc.json",
        section_name=f"Feature importance ({model_name})",
        description=description
        or "Per-feature scores for the selected reprofiler (native attributes, then permutation).",
        col1_header="Feature",
        id_field="feature",
        parent_id="samovar_feature_importance",
        parent_name="Feature importance",
        numeric_fields=("importance", "native", "permutation", "regen_agreement", "initial_agreement"),
    )
    series = {
        "Importance": {str(r["feature"]): float(r["importance"]) for r in rows},
    }
    write_bargraph_mqc(
        series,
        folder / f"feature_importance_{model_slug}_bars_mqc.json",
        section_name=f"Feature importance bars ({model_name})",
        description=description or "Relative feature importance for the selected model.",
        xlab="Importance",
        parent_id="samovar_feature_importance",
        parent_name="Feature importance",
        ymin=0,
    )
    extras = _write_optional_plots(table, folder, model_name)
    payload = {
        "model": model_name,
        "n_features": int(len(table)),
        "features": rows,
        "plots": extras,
    }
    (folder / f"feature_importance_{model_slug}.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    return payload


def score_feature_importance(
    model,
    annotation,
    initial_abundance,
    regenerated_abundance,
    config,
):
    """Contract entry: score one (or the selected) model and write MultiQC JSON."""
    cfg = dict(config or {})
    dest = _plot_dir_from_config(cfg)
    model_name = str(cfg.get("model_name") or cfg.get("best_model_name") or "best")
    names = cfg.get("feature_cols") or cfg.get("feature_names")
    table = importance_table(
        model=model,
        annotation=annotation,
        initial_abundance=initial_abundance,
        regenerated_abundance=regenerated_abundance,
        config=cfg,
        feature_names=names,
    )
    payload = write_feature_importance_outputs(
        table,
        dest,
        model_name=model_name,
        description=str(cfg.get("description") or ""),
    )
    artifact = Path(cfg.get("output_dir") or dest)
    if artifact.name == "reprofiled_annotations" or cfg.get("copy_next_to_model"):
        extra = Path(cfg.get("output_dir") or dest) / "feature_importance"
        extra.mkdir(parents=True, exist_ok=True)
        table.to_csv(extra / f"feature_importance_{_slug(model_name)}.tsv", sep="\t", index=False)
    payload["plot_dir"] = str(dest)
    return payload


def _load_imported(name: str):
    matched, spec = lookup_feature_importance(name)
    path = Path(tool_path(spec, matched)).expanduser()
    if path.suffix.lower() != ".py" or not path.is_file():
        raise MissingFeatureImportanceError(
            f"Imported feature-importance {name!r} must be a Python module."
        )
    loaded = importlib.util.spec_from_file_location(f"samovar_fi_{matched}", path)
    if loaded is None or loaded.loader is None:
        raise MissingFeatureImportanceError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(loaded)
    loaded.loader.exec_module(module)
    fn = getattr(module, "score_feature_importance", None) or getattr(module, "score", None)
    if not callable(fn):
        inst = None
        for cls_name in ("FeatureImportance", "Scorer"):
            cls = getattr(module, cls_name, None)
            if cls is not None:
                inst = cls()
                break
        if inst is not None:
            fn = getattr(inst, "score_feature_importance", None) or getattr(inst, "score", None)
    if not callable(fn):
        raise MissingFeatureImportanceError(
            f"{path} must define score_feature_importance(model, annotation, "
            "initial_abundance, regenerated_abundance, config)"
        )
    return fn


def selected_models(result: Any, *, all_models: bool) -> List[Tuple[str, Any]]:
    model = getattr(result, "model", None)
    models = getattr(result, "models", None) or {}
    if all_models and isinstance(models, dict) and models:
        out = []
        for name, est in models.items():
            if est is not None:
                out.append((str(name), est))
        if out:
            return out
    if model is not None:
        name = "best"
        if isinstance(models, dict):
            for key, est in models.items():
                if est is model:
                    name = str(key)
                    break
        return [(name, model)]
    return []


def maybe_score_reprofiler(
    result: Any,
    *,
    annotation: Any,
    initial_abundance: Any,
    regenerated_abundance: Any,
    config: Optional[Dict[str, Any]] = None,
) -> List[Dict[str, Any]]:
    cfg = dict(config or {})
    raw_name = cfg.get("feature_importance")
    if raw_name is None:
        raw_name = "builtin"
    try:
        kind = canonicalize_feature_importance(raw_name)
    except MissingFeatureImportanceError as exc:
        print(f"[feature-importance] skipped: {exc}")
        return []
    if kind == "none":
        return []
    all_models = bool(cfg.get("feature_importance_all") or cfg.get("all_models"))
    chosen = selected_models(result, all_models=all_models)
    if not chosen:
        return []
    imported_flags = ""
    if kind != "builtin":
        try:
            _matched, spec = lookup_feature_importance(kind)
            imported_flags = tool_flags(spec, kind)
        except Exception:
            imported_flags = ""
    extra = list(cfg.get("feature_importance_extra_argv") or [])
    from samovar.table_regenerators import extra_flags_argv

    extra = extra_flags_argv(cfg.get("feature_importance_flags")) + extra
    named = (cfg.get("feature_importance_tool_flags") or {}).get(kind) or (
        cfg.get("feature_importance_tool_flags") or {}
    ).get(str(raw_name))
    extra = extra + extra_flags_argv(named)
    payloads: List[Dict[str, Any]] = []
    for model_name, est in chosen:
        local = dict(cfg)
        local["model_name"] = model_name
        local["feature_cols"] = getattr(result, "feature_cols", None) or cfg.get("feature_cols")
        local["extra_argv"] = extra_flags_argv(imported_flags) + extra
        if kind == "builtin":
            payloads.append(
                score_feature_importance(
                    est, annotation, initial_abundance, regenerated_abundance, local
                )
            )
        else:
            fn = _load_imported(kind)
            payloads.append(
                fn(est, annotation, initial_abundance, regenerated_abundance, local)
            )
    try:
        out = Path(cfg.get("output_dir") or ".")
        root = out.parent if out.name == "reprofiled_annotations" else out
        from samovar.annotation_correlation import refresh_spearman_heatmaps

        refresh_spearman_heatmaps(root)
    except Exception as exc:
        print(f"[spearman] refresh skipped: {exc}", file=sys.stderr)
    return payloads


def resolve_abundance_inputs(
    config: Dict[str, Any],
    *,
    regenerated_tables: Optional[Dict[str, pd.DataFrame]] = None,
    output_dir: Optional[PathLike] = None,
) -> Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]:
    initial = _as_tables(config.get("initial_abundance"))
    if not initial:
        init_dir = config.get("initial_abundance_dir")
        if not init_dir and output_dir:
            parent = Path(output_dir)
            if parent.name == "reprofiled_annotations":
                init_dir = parent.parent / "initial_abundance"
            else:
                init_dir = parent / "initial_abundance"
        initial = _as_tables(init_dir)
    regenerated = regenerated_tables or _as_tables(config.get("regenerated_abundance"))
    if not regenerated:
        regen_dir = config.get("ground_truth_dir") or config.get("regenerated_abundance_dir")
        regenerated = _as_tables(regen_dir)
    return initial, regenerated
