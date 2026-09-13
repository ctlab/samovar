"""Sample-level abundance quality scoring.

Contract group ``sample_scoring`` (``samovar tools import --type sample-score``).

Independent of the pipeline: ``score_samples(generated, reference, config)``
returns one quality value per generated sample (higher is better).

Built-in ``bray_curtis``: for each generated sample, mean Bray–Curtis distance
to the ``k`` nearest reference samples (default k=3), then
``quality = 1 - mean_distance``.
"""

from __future__ import annotations

import argparse
import importlib.util
import math
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple, Union

import numpy as np
import pandas as pd

from samovar.abundance import (
    abundance_to_matrix,
    input_to_abundance_tables,
    is_abundance_table,
    load_abundance_dir,
    n_sample_columns,
    normalize_abundance_table,
    observed_abundance_dir,
    regenerated_abundance_dir,
)
from samovar.main_config import flags_target_matches, iter_tools, parse_tool_entry, tool_path
from samovar.paths import load_config

PathLike = Union[str, Path]
SAMPLE_SCORING_GROUP = "sample_scoring"
UNCLASSIFIED = frozenset({"0", "unclassified", "unclassified_root", "nan", "none", ""})
BUILTIN_SAMPLE_SCORERS = frozenset(
    {
        "bray_curtis",
        "bray-curtis",
        "braycurtis",
        "sample_bray",
        "nearest3_bray",
        "nearest_3_bray",
    }
)
NONE_TOKENS = frozenset({"", "none", "off", "false", "0"})
SAMPLE_QC_DIR = "sampleQC"
FULL_PHASE = "full"
FINAL_PHASE = "final"


class MissingSampleScorerError(ValueError):
    """``tools.<name>`` is missing or is not ``--type sample-score``."""


def flags_apply_to_sample_scorer(target: str, *names: Optional[str]) -> bool:
    cleaned = [name for name in names if name]
    return flags_target_matches(
        target,
        *cleaned,
        groups=(
            "sample_scoring",
            "sample-scoring",
            "sample_score",
            "sample-score",
            "sample_qc",
            "sample-qc",
        ),
    )


def canonicalize_sample_scorer(name: Optional[str]) -> str:
    key = str(name or "").strip()
    if not key or key.lower() in NONE_TOKENS:
        return ""
    low = key.lower().replace("-", "_")
    if low in {"bray_curtis", "braycurtis", "sample_bray", "nearest3_bray", "nearest_3_bray"}:
        return "bray_curtis"
    try:
        matched, _spec = lookup_sample_scorer(key)
        return matched
    except MissingSampleScorerError:
        pass
    raise ValueError(
        f"Unknown sample scorer {name!r}. Built-in: bray_curtis. "
        "Or import one with `samovar tools import --type sample-score`."
    )


def lookup_sample_scorer(name: str):
    key = str(name or "").strip()
    if not key:
        raise MissingSampleScorerError(
            "Empty sample scorer name. Import a tool with "
            "`samovar tools import -n NAME --type sample-score`."
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
        raise MissingSampleScorerError(
            f"sample-score tool {key!r} is not in the install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type sample-score`."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != SAMPLE_SCORING_GROUP:
        raise MissingSampleScorerError(
            f"tools.{matched} has group {group!r}, expected {SAMPLE_SCORING_GROUP!r}. "
            "Re-import with --type sample-score."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingSampleScorerError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return matched, parsed


def parse_sample_score_tokens(
    tokens: Optional[Sequence[str]],
) -> Tuple[str, Dict[str, str]]:
    """Split CLI/YAML tokens into (global_scorer, annotator→scorer).

    ``bray_curtis`` is global. ``kaiju:imported`` is annotator-specific.
    ``none`` clears both. Annotator keys win over the global name at resolve time.
    """
    global_name = ""
    by_annotator: Dict[str, str] = {}
    for raw in tokens or []:
        text = str(raw or "").strip()
        if not text:
            continue
        if text.lower() in NONE_TOKENS:
            global_name = ""
            by_annotator = {}
            continue
        if ":" in text:
            left, right = text.split(":", 1)
            ann, scorer = left.strip(), right.strip()
            if ann and scorer:
                if scorer.lower() in NONE_TOKENS:
                    by_annotator.pop(ann, None)
                else:
                    by_annotator[ann] = scorer
            continue
        global_name = "" if text.lower() in NONE_TOKENS else text
    return global_name, by_annotator


def ingest_sample_score_settings(
    raw: Any = None,
    by_annotator: Optional[Mapping[str, Any]] = None,
) -> Tuple[str, Dict[str, str]]:
    """Normalize YAML ``sample_score`` / ``sample_score_by_annotator`` values."""
    mapping: Dict[str, str] = {}
    tokens: List[str] = []
    if isinstance(raw, Mapping):
        mapping.update({str(k).strip(): str(v).strip() for k, v in raw.items() if str(k).strip()})
    elif isinstance(raw, (list, tuple)):
        tokens.extend(str(item) for item in raw)
    elif raw not in (None, False):
        tokens.append(str(raw))
    if isinstance(by_annotator, Mapping):
        mapping.update(
            {str(k).strip(): str(v).strip() for k, v in by_annotator.items() if str(k).strip()}
        )
    global_name, specific = parse_sample_score_tokens(tokens)
    for key, value in mapping.items():
        if not value or value.lower() in NONE_TOKENS:
            specific.pop(key, None)
        else:
            specific[key] = value
    return global_name, specific


def parse_neighbor_k_from_flags(flags: Optional[str]) -> Optional[int]:
    text = str(flags or "").strip()
    if not text:
        return None
    parts = text.replace("=", " ").split()
    for i, tok in enumerate(parts):
        if tok in {"--k", "-k", "--n-neighbors", "--n_neighbors", "--nearest"}:
            if i + 1 < len(parts):
                try:
                    return int(parts[i + 1])
                except (TypeError, ValueError):
                    return None
    return None


def sample_qc_configured(
    global_name: Optional[str] = None,
    by_annotator: Optional[Mapping[str, str]] = None,
) -> bool:
    if str(global_name or "").strip() and str(global_name).strip().lower() not in NONE_TOKENS:
        return True
    return any(str(v or "").strip() for v in (by_annotator or {}).values())


def resolve_sample_scorer(
    annotator: str,
    *,
    global_name: Optional[str] = None,
    by_annotator: Optional[Mapping[str, str]] = None,
) -> str:
    """Most specific mapping wins: annotator entry, else global, else empty."""
    specific = str((by_annotator or {}).get(annotator) or "").strip()
    if not specific:
        low = str(annotator or "").strip().lower()
        for key, value in (by_annotator or {}).items():
            if str(key).strip().lower() == low:
                specific = str(value or "").strip()
                break
    if specific:
        return canonicalize_sample_scorer(specific)
    return canonicalize_sample_scorer(global_name)


def bray_curtis_distance(left: np.ndarray, right: np.ndarray) -> float:
    a = np.asarray(left, dtype=float)
    b = np.asarray(right, dtype=float)
    a = np.where(np.isfinite(a), np.maximum(a, 0.0), 0.0)
    b = np.where(np.isfinite(b), np.maximum(b, 0.0), 0.0)
    den = float(a.sum() + b.sum())
    if den <= 0:
        return 0.0
    return float(np.abs(a - b).sum()) / den


def _drop_unclassified_rows(matrix: pd.DataFrame) -> pd.DataFrame:
    if matrix is None or matrix.empty:
        return matrix if matrix is not None else pd.DataFrame()
    keep = ~matrix.index.astype(str).str.lower().isin(UNCLASSIFIED)
    return matrix.loc[keep]


def _as_matrix(table: Any) -> pd.DataFrame:
    if table is None:
        return pd.DataFrame()
    frame: Optional[pd.DataFrame] = None
    if isinstance(table, pd.DataFrame):
        frame = table
    elif isinstance(table, (str, Path)):
        path = Path(table)
        if path.is_dir():
            loaded = load_abundance_dir(path)
            frame = next(iter(loaded.values())) if loaded else pd.DataFrame()
        elif path.is_file():
            frame = pd.read_csv(path)
        else:
            return pd.DataFrame()
    else:
        tables = input_to_abundance_tables(table)
        frame = next(iter(tables.values())) if tables else None
    if frame is None or getattr(frame, "empty", True):
        return pd.DataFrame()
    if is_abundance_table(frame) or "taxid" in getattr(frame, "columns", []) or n_sample_columns(frame):
        mat = abundance_to_matrix(normalize_abundance_table(frame))
    else:
        return pd.DataFrame()
    if mat is None or mat.empty:
        return pd.DataFrame()
    return _drop_unclassified_rows(mat)


def _align_matrices(generated: pd.DataFrame, reference: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    taxa = sorted(set(generated.index.astype(str)) | set(reference.index.astype(str)))
    gen = generated.copy()
    ref = reference.copy()
    gen.index = gen.index.astype(str)
    ref.index = ref.index.astype(str)
    gen = gen.reindex(taxa, fill_value=0.0)
    ref = ref.reindex(taxa, fill_value=0.0)
    return gen.fillna(0.0), ref.fillna(0.0)


def _neighbor_k(config: Optional[Mapping[str, Any]], n_ref: int) -> int:
    cfg = dict(config or {})
    raw = cfg.get("k", cfg.get("n_neighbors", cfg.get("nearest", 3)))
    try:
        k = int(raw)
    except (TypeError, ValueError):
        k = 3
    if k < 1:
        k = 1
    if n_ref <= 0:
        return 0
    return min(k, n_ref)


def score_samples_bray_curtis(
    generated: Any,
    reference: Any,
    config: Optional[Mapping[str, Any]] = None,
) -> pd.DataFrame:
    """Per-sample quality = 1 − mean Bray–Curtis to the k nearest reference samples."""
    gen = _as_matrix(generated)
    ref = _as_matrix(reference)
    rows: List[Dict[str, Any]] = []
    if gen.empty or not list(gen.columns):
        return pd.DataFrame(columns=["sample", "quality", "mean_distance", "n_neighbors", "scorer", "ok"])
    if ref.empty or not list(ref.columns):
        for col in gen.columns:
            rows.append(
                {
                    "sample": str(col),
                    "quality": 0.0,
                    "mean_distance": float("nan"),
                    "n_neighbors": 0,
                    "scorer": "bray_curtis",
                    "ok": False,
                }
            )
        return pd.DataFrame(rows)
    gen, ref = _align_matrices(gen, ref)
    gen_arr = gen.to_numpy(dtype=float)
    ref_arr = ref.to_numpy(dtype=float)
    k = _neighbor_k(config, ref_arr.shape[1])
    for i, col in enumerate(gen.columns):
        vec = gen_arr[:, i]
        if not np.isfinite(vec).any():
            rows.append(
                {
                    "sample": str(col),
                    "quality": float("nan"),
                    "mean_distance": float("nan"),
                    "n_neighbors": 0,
                    "scorer": "bray_curtis",
                    "ok": False,
                }
            )
            continue
        dists = np.array(
            [bray_curtis_distance(vec, ref_arr[:, j]) for j in range(ref_arr.shape[1])],
            dtype=float,
        )
        dists = dists[np.isfinite(dists)]
        if dists.size == 0 or k <= 0:
            rows.append(
                {
                    "sample": str(col),
                    "quality": 0.0,
                    "mean_distance": float("nan"),
                    "n_neighbors": 0,
                    "scorer": "bray_curtis",
                    "ok": False,
                }
            )
            continue
        nearest = np.sort(dists)[:k]
        mean = float(nearest.mean())
        quality = 1.0 - mean
        if not math.isfinite(quality):
            quality = 0.0
        quality = max(0.0, min(1.0, quality))
        rows.append(
            {
                "sample": str(col),
                "quality": quality,
                "mean_distance": mean,
                "n_neighbors": int(nearest.size),
                "scorer": "bray_curtis",
                "ok": True,
            }
        )
    return pd.DataFrame(rows)


def _normalize_score_frame(payload: Any, scorer: str) -> pd.DataFrame:
    if payload is None:
        return pd.DataFrame(columns=["sample", "quality", "scorer", "ok"])
    if isinstance(payload, pd.DataFrame):
        frame = payload.copy()
    elif isinstance(payload, Mapping):
        if "quality" in payload and "sample" in payload:
            frame = pd.DataFrame(payload)
        else:
            frame = pd.DataFrame(
                [{"sample": str(k), "quality": v} for k, v in payload.items()]
            )
    else:
        raise TypeError(f"{scorer} score_samples() must return a DataFrame or dict")
    if "sample" not in frame.columns and frame.index.name:
        frame = frame.reset_index().rename(columns={frame.columns[0]: "sample"})
    if "sample" not in frame.columns:
        raise TypeError(f"{scorer} score_samples() needs a sample column")
    if "quality" not in frame.columns:
        raise TypeError(f"{scorer} score_samples() needs a quality column")
    frame["sample"] = frame["sample"].astype(str)
    frame["quality"] = pd.to_numeric(frame["quality"], errors="coerce")
    if "scorer" not in frame.columns:
        frame["scorer"] = scorer
    if "ok" not in frame.columns:
        frame["ok"] = frame["quality"].notna()
    return frame


def _load_imported_sample_scorer(name: str) -> Any:
    matched, spec = lookup_sample_scorer(name)
    path = Path(tool_path(spec, matched))
    if path.suffix.lower() != ".py":
        raise ValueError(
            f"Imported sample scorer {name!r} must be a Python module with score_samples()."
        )
    loaded = importlib.util.spec_from_file_location(f"samovar_sample_scorer_{matched}", path)
    if loaded is None or loaded.loader is None:
        raise ValueError(f"Cannot load sample scorer {path}")
    module = importlib.util.module_from_spec(loaded)
    loaded.loader.exec_module(module)
    return module


def score_samples(
    generated: Any,
    reference: Any,
    config: Optional[Mapping[str, Any]] = None,
    *,
    scorer: Optional[str] = None,
) -> pd.DataFrame:
    """Run a sample-scoring contract. ``config`` may carry extra method fields."""
    cfg = dict(config or {})
    kind = canonicalize_sample_scorer(scorer or cfg.get("sample_score") or cfg.get("scorer") or "bray_curtis")
    if not kind:
        return pd.DataFrame(columns=["sample", "quality", "scorer", "ok"])
    if kind == "bray_curtis":
        frame = score_samples_bray_curtis(generated, reference, cfg)
        frame["scorer"] = "bray_curtis"
        return frame
    module = _load_imported_sample_scorer(kind)
    fn = getattr(module, "score_samples", None)
    if not callable(fn):
        raise TypeError(f"Imported sample scorer {kind!r} needs score_samples(generated, reference, config)")
    return _normalize_score_frame(fn(generated, reference, cfg), kind)


def sample_qc_dir(output_dir: PathLike, phase: str) -> Path:
    name = FULL_PHASE if str(phase).strip().lower() in {"full", "intermediate", "candidates"} else FINAL_PHASE
    return Path(output_dir) / SAMPLE_QC_DIR / name


def write_sample_scores(dest: PathLike, frame: pd.DataFrame, annotator: str = "") -> Path:
    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    out = frame.copy()
    if annotator and "annotator" not in out.columns:
        out.insert(0, "annotator", annotator)
    out.to_csv(path, index=False)
    return path


def _iter_generated_tables(
    output_dir: PathLike,
    phase: str,
) -> Iterable[Tuple[str, str, pd.DataFrame]]:
    """Yield ``(mode, annotator, table)`` for full (all candidates) or final (selected)."""
    from samovar.table_scorers import load_tables_by_mode_from_run

    root = Path(output_dir)
    if str(phase).strip().lower() in {FINAL_PHASE, "selected", "winner"}:
        loaded = load_abundance_dir(regenerated_abundance_dir(root))
        for annotator, table in (loaded or {}).items():
            yield FINAL_PHASE, annotator, table
        return
    tables_by_mode = load_tables_by_mode_from_run(root)
    if not tables_by_mode:
        loaded = load_abundance_dir(regenerated_abundance_dir(root))
        if loaded:
            tables_by_mode = {"generated": loaded}
    for mode, tables in (tables_by_mode or {}).items():
        for annotator, table in (tables or {}).items():
            yield str(mode), annotator, table


def stage_score_sample_qc(
    output_dir: PathLike,
    config: Optional[Mapping[str, Any]] = None,
    *,
    phase: str = FULL_PHASE,
) -> Dict[str, Any]:
    """Write per-sample scores under ``$output_dir/sampleQC/{full|final}/``."""
    cfg = dict(config or {})
    global_name = cfg.get("sample_score") or cfg.get("sample_qc") or ""
    by_annotator = dict(cfg.get("sample_score_by_annotator") or {})
    if not sample_qc_configured(global_name, by_annotator):
        return {"enabled": False, "written": []}
    root = Path(output_dir)
    observed_dir = observed_abundance_dir(root)
    observed_tables = load_abundance_dir(observed_dir)
    dest_root = sample_qc_dir(root, phase)
    dest_root.mkdir(parents=True, exist_ok=True)
    written: List[str] = []
    combined: List[pd.DataFrame] = []
    for mode, annotator, table in _iter_generated_tables(root, phase):
        scorer = resolve_sample_scorer(
            annotator, global_name=str(global_name or ""), by_annotator=by_annotator
        )
        if not scorer:
            continue
        reference = observed_tables.get(annotator)
        if reference is None and len(observed_tables) == 1:
            reference = next(iter(observed_tables.values()))
        payload = dict(cfg)
        payload["scorer"] = scorer
        payload["sample_score"] = scorer
        payload["annotator"] = annotator
        payload["phase"] = phase
        named = dict(cfg.get("sample_score_tool_flags") or {})
        flag_k = parse_neighbor_k_from_flags(
            " ".join(
                str(x)
                for x in (cfg.get("sample_score_flags"), named.get(scorer), named.get(annotator))
                if x
            )
        )
        if flag_k is not None and "k" not in payload:
            payload["k"] = flag_k
        frame = score_samples(table, reference, payload, scorer=scorer)
        frame["annotator"] = annotator
        frame["mode"] = mode
        frame["phase"] = FULL_PHASE if str(phase).lower() not in {FINAL_PHASE, "selected", "winner"} else FINAL_PHASE
        sub = dest_root / mode if mode and mode != FINAL_PHASE and str(phase).lower() not in {FINAL_PHASE, "selected", "winner"} else dest_root
        sub.mkdir(parents=True, exist_ok=True)
        path = write_sample_scores(sub / f"{annotator}.csv", frame, annotator=annotator)
        written.append(str(path))
        combined.append(frame)
    if combined:
        write_sample_scores(dest_root / "all.csv", pd.concat(combined, ignore_index=True))
        written.append(str(dest_root / "all.csv"))
    return {"enabled": True, "phase": str(phase), "written": written, "directory": str(dest_root)}


def main(argv: Optional[Sequence[str]] = None) -> int:
    import yaml

    from samovar.paths import add_output_dir_argument

    parser = argparse.ArgumentParser(prog="python -m samovar.sample_scorers")
    sub = parser.add_subparsers(dest="command", required=True)
    score = sub.add_parser("score", help="Score generated samples against a reference table")
    score.add_argument("--generated", required=True, help="Generated abundance CSV or directory")
    score.add_argument("--reference", required=True, help="Reference/predicted abundance CSV or directory")
    score.add_argument("-o", "--output", required=True, help="Output CSV")
    score.add_argument("--scorer", default="bray_curtis")
    score.add_argument("--k", type=int, default=3)
    stage = sub.add_parser("stage", help="Score a SamovaR run (sampleQC/full or sampleQC/final)")
    add_output_dir_argument(stage, required=True)
    stage.add_argument("--config", default="")
    stage.add_argument("--phase", default=FULL_PHASE, choices=[FULL_PHASE, FINAL_PHASE])
    args = parser.parse_args(list(argv) if argv is not None else None)
    if args.command == "score":
        frame = score_samples(
            args.generated,
            args.reference,
            {"k": args.k, "scorer": args.scorer},
            scorer=args.scorer,
        )
        write_sample_scores(args.output, frame)
        print(f"wrote {args.output} n={len(frame)}")
        return 0
    cfg: Dict[str, Any] = {}
    if args.config:
        cfg = yaml.safe_load(Path(args.config).read_text(encoding="utf-8")) or {}
    result = stage_score_sample_qc(args.output_dir, cfg, phase=args.phase)
    print(f"sampleQC phase={result.get('phase')} written={len(result.get('written') or [])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
