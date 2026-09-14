"""Sample-level abundance filtering (uses sample quality scores).

Contract group ``sample_filtering`` (``samovar tools import --type sample-filter``).

Independent of the pipeline: ``filter_samples(table, scores, config)`` drops
generated samples. Higher quality is better. Built-ins:

* ``top``: keep top N and/or top fraction of samples.
* ``median_sd``: keep samples with quality in median ± X·SD.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
import warnings
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple, Union

import pandas as pd

from samovar.abundance import (
    load_abundance_dir,
    n_sample_columns,
    observed_abundance_dir,
    regenerated_abundance_dir,
    write_abundance_dir,
)
from samovar.main_config import flags_target_matches, iter_tools, parse_tool_entry, tool_path
from samovar.paths import load_config
from samovar.sample_scorers import (
    NONE_TOKENS,
    FINAL_PHASE,
    FULL_PHASE,
    _iter_generated_tables,
    canonicalize_sample_scorer,
    resolve_sample_scorer,
    sample_qc_configured,
    sample_qc_dir,
    score_samples,
    stage_score_sample_qc,
)
from samovar.table_regenerators import is_table_method_name

PathLike = Union[str, Path]
SAMPLE_FILTERING_GROUP = "sample_filtering"
SAMPLE_FILTER_DIR = "sampleFilter"
KEEP_ONE_WARNING = (
    "filtration remove all samples; keeping the 1 best regenerated sample"
)

BUILTIN_FILTERS = {
    "top": "top",
    "top_n": "top",
    "keep_top": "top",
    "topk": "top",
    "median_sd": "median_sd",
    "median-sd": "median_sd",
    "sd": "median_sd",
    "sigma": "median_sd",
    "median": "median_sd",
}

FILTER_PARAM_KEYS = {
    "n": ("n", "top_n", "keep_n"),
    "frac": ("frac", "fraction", "percent", "p", "pct"),
    "sd": ("sd", "n_sd", "sigma", "x", "nsd"),
}

CLI_PARAM_FLAGS = {
    "n": "--sample-filter-n",
    "frac": "--sample-filter-frac",
    "sd": "--sample-filter-sd",
}

BUILTIN_SUPPORTED_PARAMS = {
    "top": frozenset({"n", "frac"}),
    "median_sd": frozenset({"sd"}),
}

STAGE_CANDIDATES = frozenset(
    {"candidates", "full", "before", "before_compare", "intermediate"}
)
STAGE_FINAL = frozenset({"final", "selected", "winner"})
STAGE_BOTH = frozenset({"both", "all"})


class MissingSampleFilterError(ValueError):
    """``tools.<name>`` is missing or is not ``--type sample-filter``."""


def flags_apply_to_sample_filter(target: str, *names: Optional[str]) -> bool:
    cleaned = [name for name in names if name]
    return flags_target_matches(
        target,
        *cleaned,
        groups=(
            "sample_filtering",
            "sample-filtering",
            "sample_filter",
            "sample-filter",
            "sample_qc_filter",
        ),
    )


def canonicalize_sample_filter(name: Optional[str]) -> str:
    key = str(name or "").strip()
    if not key or key.lower() in NONE_TOKENS:
        return ""
    low = key.lower().replace("-", "_")
    if low in BUILTIN_FILTERS:
        return BUILTIN_FILTERS[low]
    try:
        matched, _spec = lookup_sample_filter(key)
        return matched
    except MissingSampleFilterError:
        pass
    raise ValueError(
        f"Unknown sample filter {name!r}. Built-in: top, median_sd. "
        "Or import one with `samovar tools import --type sample-filter`."
    )


def lookup_sample_filter(name: str):
    key = str(name or "").strip()
    if not key:
        raise MissingSampleFilterError(
            "Empty sample filter name. Import a tool with "
            "`samovar tools import -n NAME --type sample-filter`."
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
        raise MissingSampleFilterError(
            f"sample-filter tool {key!r} is not in the install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type sample-filter`."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != SAMPLE_FILTERING_GROUP:
        raise MissingSampleFilterError(
            f"tools.{matched} has group {group!r}, expected {SAMPLE_FILTERING_GROUP!r}. "
            "Re-import with --type sample-filter."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingSampleFilterError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return matched, parsed


def parse_sample_filter_stage(raw: Any) -> str:
    text = str(raw or "final").strip().lower().replace("-", "_")
    if text in STAGE_BOTH:
        return "both"
    if text in STAGE_CANDIDATES:
        return "candidates"
    if text in STAGE_FINAL or not text:
        return "final"
    raise ValueError(
        f"Unknown sample-filter stage {raw!r}. Use final, candidates, or both."
    )


def sample_filter_configured(
    global_name: Optional[str] = None,
    by_method: Optional[Mapping[str, str]] = None,
) -> bool:
    if str(global_name or "").strip() and str(global_name).strip().lower() not in NONE_TOKENS:
        return True
    return any(str(v or "").strip() for v in (by_method or {}).values())


def filter_candidates_enabled(stage: str) -> bool:
    return parse_sample_filter_stage(stage) in {"candidates", "both"}


def filter_final_enabled(stage: str) -> bool:
    return parse_sample_filter_stage(stage) in {"final", "both"}


def parse_named_tokens(
    tokens: Optional[Sequence[str]],
    *,
    table_methods: Optional[Sequence[str]] = None,
) -> Tuple[str, Dict[str, str]]:
    """Split ``NAME`` / ``method:NAME`` tokens. Methods are table regenerators."""
    global_name = ""
    by_method: Dict[str, str] = {}
    for raw in tokens or []:
        text = str(raw or "").strip()
        if not text:
            continue
        if text.lower() in NONE_TOKENS:
            global_name = ""
            by_method = {}
            continue
        if ":" in text:
            left, right = text.split(":", 1)
            key, value = left.strip(), right.strip()
            if key and value:
                if not is_table_method_name(key, list(table_methods or [])):
                    raise ValueError(
                        f"Unknown table generation method {key!r} in sample-filter. "
                        "Use a regeneration method (direct, bootstrap, vae, glm, …)."
                    )
                if value.lower() in NONE_TOKENS:
                    by_method.pop(key, None)
                else:
                    by_method[key] = value
            continue
        global_name = "" if text.lower() in NONE_TOKENS else text
    # Drop annotator-shaped keys that are not table methods when extras known.
    if table_methods is not None:
        cleaned: Dict[str, str] = {}
        for key, value in by_method.items():
            if is_table_method_name(key, list(table_methods)):
                cleaned[key] = value
            else:
                cleaned[key] = value
        by_method = cleaned
    return global_name, by_method


def parse_param_tokens(
    tokens: Optional[Sequence[Any]],
) -> Tuple[Optional[float], Dict[str, float]]:
    """``5`` is global; ``bootstrap:8`` is per-method."""
    global_value: Optional[float] = None
    by_method: Dict[str, float] = {}
    for raw in tokens or []:
        text = str(raw or "").strip()
        if not text or text.lower() in NONE_TOKENS:
            global_value = None
            by_method = {}
            continue
        if ":" in text:
            left, right = text.split(":", 1)
            key, value = left.strip(), right.strip()
            if key and value:
                by_method[key] = float(value)
            continue
        global_value = float(text)
    return global_value, by_method


def _lookup_ci(mapping: Optional[Mapping[str, Any]], key: str) -> Any:
    if not mapping:
        return None
    if key in mapping:
        return mapping.get(key)
    low = str(key).strip().lower()
    for stored, value in mapping.items():
        if str(stored).strip().lower() == low:
            return value
    return None


def resolve_sample_filter(
    method: str,
    *,
    global_name: Optional[str] = None,
    by_method: Optional[Mapping[str, str]] = None,
) -> str:
    specific = str(_lookup_ci(by_method, method) or "").strip()
    if specific:
        return canonicalize_sample_filter(specific)
    return canonicalize_sample_filter(global_name)


def _param_from_config(config: Optional[Mapping[str, Any]], kind: str) -> Optional[float]:
    cfg = dict(config or {})
    for key in FILTER_PARAM_KEYS[kind]:
        if cfg.get(key) not in (None, ""):
            try:
                return float(cfg.get(key))
            except (TypeError, ValueError):
                raise ValueError(f"Invalid sample filter {kind}={cfg.get(key)!r}")
    return None


def collect_filter_params(config: Optional[Mapping[str, Any]]) -> Dict[str, float]:
    out: Dict[str, float] = {}
    for kind in ("n", "frac", "sd"):
        value = _param_from_config(config, kind)
        if value is not None:
            out[kind] = value
    return out


def supported_params_for(filter_name: str) -> frozenset:
    kind = canonicalize_sample_filter(filter_name)
    if kind in BUILTIN_SUPPORTED_PARAMS:
        return BUILTIN_SUPPORTED_PARAMS[kind]
    try:
        module = _load_imported_sample_filter(kind)
    except Exception:
        return frozenset()
    raw = getattr(module, "SUPPORTED_PARAMS", None) or getattr(module, "supported_params", None)
    if callable(raw):
        raw = raw()
    if not raw:
        return frozenset()
    return frozenset(str(x).strip().lower() for x in raw if str(x).strip())


def _supported_cli_flags(params: frozenset) -> str:
    flags = [CLI_PARAM_FLAGS[k] for k in ("n", "frac", "sd") if k in params]
    if not flags:
        return "a supported parameter"
    if len(flags) == 1:
        return flags[0]
    return " or ".join(flags)


def unsupported_param_message(filter_name: str, bad: Sequence[str]) -> str:
    supported = supported_params_for(filter_name)
    hint = _supported_cli_flags(supported)
    bad_flags = " / ".join(CLI_PARAM_FLAGS.get(k, k) for k in bad)
    _ = bad_flags
    return (
        "Not supported by the exact abundance table filtering function. "
        f"Use other {hint} or filtering function"
    )


def validate_filter_params_at_prepare(
    filter_name: str,
    params: Mapping[str, Any],
    *,
    n_generated: Optional[int] = None,
) -> None:
    """Raise ValueError during ``samovar prepare`` for illegal combinations."""
    kind = canonicalize_sample_filter(filter_name)
    if not kind:
        return
    supported = supported_params_for(kind)
    present = {k: v for k, v in dict(params or {}).items() if v not in (None, "")}
    bad = [k for k in present if k not in supported]
    if bad:
        raise ValueError(unsupported_param_message(kind, bad))
    if "n" in present:
        try:
            n_val = int(present["n"])
        except (TypeError, ValueError):
            raise ValueError("top N need to keep more than 1 sample")
        if n_val <= 1:
            raise ValueError("top N need to keep more than 1 sample")
        if n_generated is not None and n_val >= int(n_generated):
            raise ValueError(
                "top N need to be smaller than the amount of the generated samples"
            )
    if "frac" in present:
        frac = float(present["frac"])
        if not (0.0 < frac < 1.0):
            raise ValueError("top % need to be more than 0 and less than 1")
    if "sd" in present:
        sd = float(present["sd"])
        if sd < 0:
            raise ValueError("sample-filter-sd must be >= 0")


def validate_pipeline_sample_filters(
    *,
    sample_filter: str = "",
    sample_filter_by_method: Optional[Mapping[str, str]] = None,
    params: Optional[Mapping[str, Any]] = None,
    params_by_method: Optional[Mapping[str, Mapping[str, Any]]] = None,
    table_methods: Optional[Sequence[str]] = None,
    n_generated: Optional[int] = None,
    sample_score: str = "",
    sample_score_by_annotator: Optional[Mapping[str, str]] = None,
    sample_score_by_method: Optional[Mapping[str, str]] = None,
) -> None:
    if not sample_filter_configured(sample_filter, sample_filter_by_method):
        return
    if not sample_qc_configured(
        sample_score, sample_score_by_annotator, sample_score_by_method
    ):
        raise ValueError(
            "Sample filtering requires a regenerated sample quality function. "
            "Pass --sample-score."
        )
    global_params = {k: v for k, v in dict(params or {}).items() if v not in (None, "")}
    per = dict(params_by_method or {})
    methods = list(table_methods or [])
    if sample_filter_by_method:
        methods = list(dict.fromkeys([*methods, *sample_filter_by_method.keys()]))
    if not methods:
        methods = [""]
    seen = set()
    for method in methods:
        filt = resolve_sample_filter(
            method or "",
            global_name=sample_filter,
            by_method=sample_filter_by_method,
        )
        if not filt:
            continue
        extra = dict(_lookup_ci(per, method) or {}) if method else {}
        merged = {**global_params, **extra}
        key = (filt, tuple(sorted(merged.items())))
        if key in seen:
            continue
        seen.add(key)
        validate_filter_params_at_prepare(filt, merged, n_generated=n_generated)
    if sample_filter and not seen:
        validate_filter_params_at_prepare(
            sample_filter, global_params, n_generated=n_generated
        )


def _scores_frame(scores: Any) -> pd.DataFrame:
    if scores is None:
        return pd.DataFrame(columns=["sample", "quality"])
    if isinstance(scores, pd.DataFrame):
        frame = scores.copy()
    elif isinstance(scores, Mapping):
        frame = pd.DataFrame(
            [{"sample": str(k), "quality": v} for k, v in scores.items()]
        )
    else:
        raise TypeError("scores must be a DataFrame or dict of sample → quality")
    if "sample" not in frame.columns:
        frame = frame.reset_index().rename(columns={frame.columns[0]: "sample"})
    if "quality" not in frame.columns:
        raise TypeError("scores need a quality column")
    frame["sample"] = frame["sample"].astype(str)
    frame["quality"] = pd.to_numeric(frame["quality"], errors="coerce")
    return frame


def sample_column_key(col: Any) -> str:
    text = str(col)
    if text.startswith("N_"):
        return text[2:]
    return text


def subset_abundance_samples(table: pd.DataFrame, keep: Sequence[str]) -> pd.DataFrame:
    keep_set = {str(x) for x in keep}
    cols: List[Any] = []
    for col in table.columns:
        name = str(col)
        if name.lower() == "taxid" or name == "taxid":
            cols.append(col)
            continue
        if sample_column_key(col) in keep_set or name in keep_set:
            cols.append(col)
    if not cols:
        return table.iloc[:, 0:0]
    return table.loc[:, cols].copy()


def _ranked_samples(scores: pd.DataFrame) -> pd.DataFrame:
    frame = scores.dropna(subset=["quality"]).copy()
    if frame.empty:
        frame = scores.copy()
        frame["quality"] = frame["quality"].fillna(float("-inf"))
    return frame.sort_values("quality", ascending=False, kind="mergesort")


def _keep_best_one(ranked: pd.DataFrame, warnings_out: List[str]) -> List[str]:
    warnings_out.append(KEEP_ONE_WARNING)
    warnings.warn(KEEP_ONE_WARNING, UserWarning, stacklevel=3)
    if ranked.empty:
        return []
    return [str(ranked.iloc[0]["sample"])]


def filter_top(
    table: pd.DataFrame,
    scores: pd.DataFrame,
    config: Optional[Mapping[str, Any]] = None,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    params = collect_filter_params(config)
    ranked = _ranked_samples(scores)
    n_samples = len(ranked)
    notes: List[str] = []
    if n_samples <= 0:
        return subset_abundance_samples(table, []), [], notes
    keep_n = n_samples
    if "n" in params:
        n_val = int(params["n"])
        if n_val <= 1:
            raise ValueError("top N need to keep more than 1 sample")
        if n_val >= n_samples:
            raise ValueError(
                "top N need to be smaller than the amount of the generated samples"
            )
        keep_n = min(keep_n, n_val)
    if "frac" in params:
        frac = float(params["frac"])
        if not (0.0 < frac < 1.0):
            raise ValueError("top % need to be more than 0 and less than 1")
        keep_n = min(keep_n, max(0, int(math.floor(frac * n_samples))))
    if keep_n <= 1:
        keep = _keep_best_one(ranked, notes)
        return subset_abundance_samples(table, keep), keep, notes
    keep = [str(s) for s in ranked["sample"].head(int(keep_n)).tolist()]
    return subset_abundance_samples(table, keep), keep, notes


def filter_median_sd(
    table: pd.DataFrame,
    scores: pd.DataFrame,
    config: Optional[Mapping[str, Any]] = None,
) -> Tuple[pd.DataFrame, List[str], List[str]]:
    params = collect_filter_params(config)
    if "sd" not in params:
        raise ValueError(unsupported_param_message("median_sd", ["n", "frac"]))
    x = float(params["sd"])
    if x < 0:
        raise ValueError("sample-filter-sd must be >= 0")
    ranked = _ranked_samples(scores)
    notes: List[str] = []
    values = pd.to_numeric(ranked["quality"], errors="coerce").dropna()
    if values.empty:
        keep = _keep_best_one(ranked, notes)
        return subset_abundance_samples(table, keep), keep, notes
    median = float(values.median())
    std = float(values.std(ddof=1)) if len(values) > 1 else 0.0
    if not math.isfinite(std) or std == 0.0:
        lo, hi = median, median
    else:
        lo, hi = median - x * std, median + x * std
    mask = ranked["quality"].between(lo, hi, inclusive="both")
    keep = [str(s) for s in ranked.loc[mask, "sample"].tolist()]
    if not keep:
        keep = _keep_best_one(ranked, notes)
    elif len(keep) == 1 and len(ranked) > 1:
        notes.append(KEEP_ONE_WARNING)
        warnings.warn(KEEP_ONE_WARNING, UserWarning, stacklevel=3)
    return subset_abundance_samples(table, keep), keep, notes


def _load_imported_sample_filter(name: str) -> Any:
    matched, spec = lookup_sample_filter(name)
    path = Path(tool_path(spec, matched))
    if path.suffix.lower() != ".py":
        raise ValueError(
            f"Imported sample filter {name!r} must be a Python module with filter_samples()."
        )
    loaded = importlib.util.spec_from_file_location(f"samovar_sample_filter_{matched}", path)
    if loaded is None or loaded.loader is None:
        raise ValueError(f"Cannot load sample filter {path}")
    module = importlib.util.module_from_spec(loaded)
    loaded.loader.exec_module(module)
    return module


def filter_samples(
    table: Any,
    scores: Any,
    config: Optional[Mapping[str, Any]] = None,
    *,
    filt: Optional[str] = None,
) -> pd.DataFrame:
    """Run a sample-filtering contract. Returns a subset abundance table."""
    cfg = dict(config or {})
    kind = canonicalize_sample_filter(
        filt or cfg.get("sample_filter") or cfg.get("filter") or ""
    )
    if isinstance(table, (str, Path)):
        table = pd.read_csv(table)
    scores_frame = _scores_frame(scores)
    if not kind:
        return table
    if kind == "top":
        out, _keep, _notes = filter_top(table, scores_frame, cfg)
        return out
    if kind == "median_sd":
        out, _keep, _notes = filter_median_sd(table, scores_frame, cfg)
        return out
    module = _load_imported_sample_filter(kind)
    fn = getattr(module, "filter_samples", None)
    if not callable(fn):
        raise TypeError(
            f"Imported sample filter {kind!r} needs filter_samples(table, scores, config)"
        )
    result = fn(table, scores_frame, cfg)
    if isinstance(result, pd.DataFrame):
        return result
    if isinstance(result, Mapping) and "table" in result:
        return result["table"]
    raise TypeError(f"{kind} filter_samples() must return a DataFrame")


def filter_samples_with_meta(
    table: pd.DataFrame,
    scores: Any,
    config: Optional[Mapping[str, Any]] = None,
    *,
    filt: Optional[str] = None,
) -> Dict[str, Any]:
    cfg = dict(config or {})
    kind = canonicalize_sample_filter(
        filt or cfg.get("sample_filter") or cfg.get("filter") or ""
    )
    scores_frame = _scores_frame(scores)
    notes: List[str] = []
    keep: List[str] = []
    if kind == "top":
        out, keep, notes = filter_top(table, scores_frame, cfg)
    elif kind == "median_sd":
        out, keep, notes = filter_median_sd(table, scores_frame, cfg)
    else:
        out = filter_samples(table, scores_frame, cfg, filt=kind)
        keep = [sample_column_key(c) for c in n_sample_columns(out)]
    dropped = [
        str(s)
        for s in scores_frame["sample"].astype(str).tolist()
        if str(s) not in set(keep)
    ]
    return {"table": out, "kept": keep, "dropped": dropped, "notes": notes, "filter": kind}


def _method_params(
    cfg: Mapping[str, Any],
    method: str,
) -> Dict[str, Any]:
    params = {
        k: cfg[k]
        for k in ("n", "top_n", "keep_n", "frac", "fraction", "percent", "sd", "n_sd", "sigma", "x")
        if cfg.get(k) not in (None, "")
    }
    by = dict(cfg.get("sample_filter_params_by_method") or {})
    extra = _lookup_ci(by, method) or {}
    params.update(dict(extra))
    # per-method scalar maps
    for kind, cfg_key in (("n", "sample_filter_n_by_method"), ("frac", "sample_filter_frac_by_method"), ("sd", "sample_filter_sd_by_method")):
        named = dict(cfg.get(cfg_key) or {})
        hit = _lookup_ci(named, method)
        if hit not in (None, ""):
            params[kind] = hit
    return params


def _load_scores_for(
    dest_root: Path,
    mode: str,
    annotator: str,
    table: pd.DataFrame,
    reference: Any,
    scorer: str,
    cfg: Mapping[str, Any],
) -> pd.DataFrame:
    candidates = [
        dest_root / mode / f"{annotator}.csv",
        dest_root / f"{annotator}.csv",
    ]
    for path in candidates:
        if path.is_file():
            return pd.read_csv(path)
    payload = dict(cfg)
    payload["scorer"] = scorer
    return score_samples(table, reference, payload, scorer=scorer)


def stage_filter_sample_qc(
    output_dir: PathLike,
    config: Optional[Mapping[str, Any]] = None,
    *,
    phase: str = FULL_PHASE,
) -> Dict[str, Any]:
    """Filter generated abundance tables using sample quality scores."""
    cfg = dict(config or {})
    global_filt = cfg.get("sample_filter") or cfg.get("sample_qc_filter") or ""
    by_method = dict(cfg.get("sample_filter_by_method") or {})
    if not sample_filter_configured(global_filt, by_method):
        return {"enabled": False, "written": []}
    stage = parse_sample_filter_stage(cfg.get("sample_filter_stage") or "final")
    want_candidates = filter_candidates_enabled(stage)
    want_final = filter_final_enabled(stage)
    phase_key = FULL_PHASE if str(phase).strip().lower() in STAGE_CANDIDATES | {FULL_PHASE} else FINAL_PHASE
    if phase_key == FULL_PHASE and not want_candidates:
        return {"enabled": False, "skipped": "stage", "written": []}
    if phase_key == FINAL_PHASE and not want_final:
        return {"enabled": False, "skipped": "stage", "written": []}
    root = Path(output_dir)
    if phase_key == FULL_PHASE:
        stage_score_sample_qc(root, cfg, phase=FULL_PHASE)
    else:
        stage_score_sample_qc(root, cfg, phase=FINAL_PHASE)
    observed = load_abundance_dir(observed_abundance_dir(root))
    score_root = sample_qc_dir(root, phase_key)
    dest_root = Path(root) / SAMPLE_FILTER_DIR / (
        FULL_PHASE if phase_key == FULL_PHASE else FINAL_PHASE
    )
    dest_root.mkdir(parents=True, exist_ok=True)
    written: List[str] = []
    reports: List[Dict[str, Any]] = []
    grouped: Dict[str, Dict[str, pd.DataFrame]] = {}
    for mode, annotator, table in _iter_generated_tables(root, phase_key):
        filt = resolve_sample_filter(
            mode, global_name=str(global_filt or ""), by_method=by_method
        )
        if not filt:
            grouped.setdefault(mode, {})[annotator] = table
            continue
        scorer = resolve_sample_scorer(
            annotator,
            method=mode,
            global_name=str(cfg.get("sample_score") or ""),
            by_annotator=dict(cfg.get("sample_score_by_annotator") or {}),
            by_method=dict(cfg.get("sample_score_by_method") or {}),
        )
        reference = observed.get(annotator)
        if reference is None and len(observed) == 1:
            reference = next(iter(observed.values()))
        scores = _load_scores_for(score_root, mode, annotator, table, reference, scorer, cfg)
        payload = dict(cfg)
        payload.update(_method_params(cfg, mode))
        payload["sample_filter"] = filt
        payload["filter"] = filt
        meta = filter_samples_with_meta(table, scores, payload, filt=filt)
        grouped.setdefault(mode, {})[annotator] = meta["table"]
        report = {
            "mode": mode,
            "annotator": annotator,
            "filter": filt,
            "kept": meta["kept"],
            "dropped": meta["dropped"],
            "notes": meta["notes"],
            "n_kept": len(meta["kept"]),
            "n_dropped": len(meta["dropped"]),
        }
        reports.append(report)
        sub = dest_root / mode if phase_key == FULL_PHASE and mode not in {FINAL_PHASE} else dest_root
        sub.mkdir(parents=True, exist_ok=True)
        path = sub / f"{annotator}.csv"
        meta["table"].to_csv(path, index=False)
        written.append(str(path))
        pd.DataFrame(
            {
                "sample": meta["kept"] + meta["dropped"],
                "kept": [True] * len(meta["kept"]) + [False] * len(meta["dropped"]),
            }
        ).to_csv(sub / f"{annotator}.kept.csv", index=False)
    if phase_key == FULL_PHASE:
        dest = regenerated_abundance_dir(root)
        candidates_root = dest / ".table_candidates"
        multi = len(grouped) > 1 or candidates_root.is_dir()
        for mode, tables in grouped.items():
            if multi:
                write_abundance_dir(candidates_root / mode, tables)
            else:
                write_abundance_dir(dest, tables)
    else:
        mixed: Dict[str, pd.DataFrame] = {}
        for tables in grouped.values():
            mixed.update(tables)
        if mixed:
            write_abundance_dir(regenerated_abundance_dir(root), mixed)
    if reports:
        pd.DataFrame(reports).to_csv(dest_root / "report.csv", index=False)
        written.append(str(dest_root / "report.csv"))
    return {
        "enabled": True,
        "phase": phase_key,
        "written": written,
        "directory": str(dest_root),
        "reports": reports,
    }


def ingest_sample_filter_settings(
    raw: Any = None,
    by_method: Optional[Mapping[str, Any]] = None,
    *,
    table_methods: Optional[Sequence[str]] = None,
) -> Tuple[str, Dict[str, str]]:
    mapping: Dict[str, str] = {}
    tokens: List[str] = []
    if isinstance(raw, Mapping):
        mapping.update({str(k).strip(): str(v).strip() for k, v in raw.items() if str(k).strip()})
    elif isinstance(raw, (list, tuple)):
        tokens.extend(str(item) for item in raw)
    elif raw not in (None, False):
        tokens.append(str(raw))
    if isinstance(by_method, Mapping):
        mapping.update(
            {str(k).strip(): str(v).strip() for k, v in by_method.items() if str(k).strip()}
        )
    global_name, specific = parse_named_tokens(tokens, table_methods=table_methods)
    for key, value in mapping.items():
        if not value or value.lower() in NONE_TOKENS:
            specific.pop(key, None)
        else:
            specific[key] = value
    return global_name, specific


def main(argv: Optional[Sequence[str]] = None) -> int:
    import yaml

    from samovar.paths import add_output_dir_argument

    parser = argparse.ArgumentParser(prog="python -m samovar.sample_filters")
    sub = parser.add_subparsers(dest="command", required=True)
    filt = sub.add_parser("filter", help="Filter one abundance table using sample scores")
    filt.add_argument("--table", required=True)
    filt.add_argument("--scores", required=True)
    filt.add_argument("-o", "--output", required=True)
    filt.add_argument("--filter", dest="filt", default="top")
    filt.add_argument("--n", type=float, default=None)
    filt.add_argument("--frac", type=float, default=None)
    filt.add_argument("--sd", type=float, default=None)
    stage = sub.add_parser("stage", help="Filter a SamovaR run (sampleFilter/full or final)")
    add_output_dir_argument(stage, required=True)
    stage.add_argument("--config", default="")
    stage.add_argument("--phase", default=FULL_PHASE, choices=[FULL_PHASE, FINAL_PHASE])
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.command == "filter":
        cfg: Dict[str, Any] = {"sample_filter": args.filt}
        if args.n is not None:
            cfg["n"] = args.n
        if args.frac is not None:
            cfg["frac"] = args.frac
        if args.sd is not None:
            cfg["sd"] = args.sd
        table = pd.read_csv(args.table)
        scores = pd.read_csv(args.scores)
        out = filter_samples(table, scores, cfg, filt=args.filt)
        Path(args.output).parent.mkdir(parents=True, exist_ok=True)
        out.to_csv(args.output, index=False)
        print(f"wrote {args.output} samples={len(n_sample_columns(out))}")
        return 0
    cfg = {}
    if args.config:
        cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
    result = stage_filter_sample_qc(args.output_dir, cfg, phase=args.phase)
    print(f"sampleFilter phase={result.get('phase')} written={len(result.get('written') or [])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
