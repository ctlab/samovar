"""Custom scoring / visualization wrappers (``samovar tools import --type scoring``).

Built-in OPAL and MultiQC keep their native CLIs and pipeline steps. Imported
tools with group ``scoring`` run after each viz checkpoint. The 6th ``tools.*``
slot is an input glob under the run output directory (default ``*annotations``).
"""

from __future__ import annotations

import argparse
import glob as globmod
import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Union

import yaml

from samovar.main_config import (
    DEFAULT_SCORING_INPUTS,
    flags_target_matches,
    iter_tools,
    merge_flag_strings,
    parse_tool_entry,
    tool_flags,
    tool_inputs,
    tool_path,
)
from samovar.paths import add_output_dir_argument, load_config
from samovar.table_regenerators import extra_flags_argv

SCORING_GROUP = "scoring"
BUILTIN_SCORING_NAMES = {"opal", "opal.py", "multiqc"}
SCORING_FLAG_GROUPS = (
    "scoring",
    "score",
    "viz",
    "visualisation",
    "visualization",
    "visualizations",
    "visualisations",
    "plots",
)

INPUT_ALIASES = {
    "annotations": "*annotations",
    "*annotation": "*annotations",
    "*_annotation": "*annotations",
    "annotation_table": "*annotations/combined_annotation_table.csv",
    "*annotation_table": "*annotations/combined_annotation_table.csv",
    "combined_annotation_table.csv": "*annotations/combined_annotation_table.csv",
    "plots": "*_annotations_plots",
    "*plots": "*_plots",
    "*_plots": "*_plots",
}

PathLike = Union[str, Path]


class MissingScorerError(ValueError):
    """``tools.<name>`` is missing or is not a scoring tool."""


def canonicalize_inputs_glob(pattern: Optional[str]) -> str:
    raw = str(pattern or "").strip()
    if not raw:
        return DEFAULT_SCORING_INPUTS
    mapped = INPUT_ALIASES.get(raw) or INPUT_ALIASES.get(raw.lower())
    return mapped or raw


def split_input_patterns(pattern: Optional[str]) -> List[str]:
    text = canonicalize_inputs_glob(pattern)
    parts: List[str] = []
    for chunk in text.replace(";", ",").split(","):
        piece = chunk.strip()
        if piece:
            parts.append(canonicalize_inputs_glob(piece))
    return parts or [DEFAULT_SCORING_INPUTS]


def expand_scoring_inputs(
    run_dir: PathLike,
    pattern: Optional[str] = None,
) -> List[Path]:
    """Resolve comma-separated globs against the run output directory."""
    root = Path(run_dir)
    found: List[Path] = []
    seen = set()
    for pat in split_input_patterns(pattern):
        hits = _glob_under(root, pat)
        if not hits and "*" not in pat and "?" not in pat and "[" not in pat:
            direct = (root / pat).expanduser()
            if direct.exists():
                hits = [direct]
        for hit in hits:
            key = str(hit.resolve()) if hit.exists() else str(hit)
            if key in seen:
                continue
            seen.add(key)
            found.append(hit)
    return found


def _glob_under(root: Path, pat: str) -> List[Path]:
    pat = pat.lstrip("/")
    query = str(root / pat)
    matches = [Path(p) for p in sorted(globmod.glob(query, recursive=True))]
    if matches:
        return matches
    try:
        return sorted(p for p in root.glob(pat) if p.exists())
    except OSError:
        return []


def flags_apply_to_scorer(target: str, name: Optional[str] = None) -> bool:
    names = [name] if name else []
    return flags_target_matches(target, *names, groups=SCORING_FLAG_GROUPS)


def is_scoring_flag_target(target: str, names: Optional[Sequence[str]] = None) -> bool:
    extra = list(names or [])
    extra.extend(iter_custom_scoring_names())
    return flags_target_matches(target, *extra, groups=SCORING_FLAG_GROUPS)


def iter_custom_scoring_names(tools: Optional[Dict[str, List[str]]] = None) -> List[str]:
    mapping = tools if tools is not None else iter_tools(load_config())
    names: List[str] = []
    for name, spec in mapping.items():
        parsed = parse_tool_entry(spec, name)
        group = str(parsed[3] or "").strip()
        if group != SCORING_GROUP:
            continue
        if name.lower() in BUILTIN_SCORING_NAMES:
            continue
        if not tool_path(parsed, name):
            continue
        names.append(name)
    return names


def lookup_scorer(name: str) -> list:
    key = str(name or "").strip()
    if not key:
        raise MissingScorerError(
            "Empty scoring tool name. Import with "
            "`samovar tools import -n NAME --type scoring --inputs '*annotations'`."
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
        raise MissingScorerError(
            f"scoring tool {key!r} is not in the install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type scoring`."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != SCORING_GROUP:
        raise MissingScorerError(
            f"tools.{matched} has group {group!r}, expected {SCORING_GROUP!r}. "
            "Re-import with --type scoring (or --type viz)."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingScorerError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return parsed


def attach_scorer_flags(name: str, spec: list, config: Dict[str, Any]) -> Dict[str, Any]:
    cfg = dict(config or {})
    named_map = cfg.get("scoring_tool_flags") or {}
    named = named_map.get(name)
    if named is None:
        low = str(name).lower()
        for key, value in named_map.items():
            if str(key).lower() == low:
                named = value
                break
    cfg["extra_flags"] = merge_flag_strings(
        tool_flags(spec, name),
        cfg.get("scoring_flags"),
        named,
        cfg.get("extra_flags"),
    )
    cfg["extra_argv"] = extra_flags_argv(cfg["extra_flags"])
    return cfg


def run_scorer(
    name: str,
    *,
    run_dir: PathLike,
    config: Optional[Dict[str, Any]] = None,
    inputs: Optional[Sequence[PathLike]] = None,
) -> Path:
    spec = lookup_scorer(name)
    cfg = attach_scorer_flags(name, spec, dict(config or {}))
    root = Path(run_dir)
    pattern = tool_inputs(spec, name)
    cfg["inputs_glob"] = pattern
    cfg["name"] = name
    resolved = [Path(p) for p in inputs] if inputs is not None else expand_scoring_inputs(
        root, pattern
    )
    dest = Path(cfg.get("score_output") or (root / f"{name}_scores"))
    dest.mkdir(parents=True, exist_ok=True)
    cfg["output_dir"] = str(dest)
    cfg["run_dir"] = str(root)
    path = Path(tool_path(spec, name)).expanduser()
    py_fn = _try_python_callable(path, name)
    if py_fn is not None:
        py_fn(resolved, dest, cfg)
        return dest
    _run_cli_scorer(path, resolved, dest, cfg)
    return dest


def run_custom_scorers(
    run_dir: PathLike,
    *,
    config: Optional[Dict[str, Any]] = None,
    names: Optional[Sequence[str]] = None,
    stage: str = "",
) -> List[Path]:
    cfg = dict(config or {})
    if stage:
        cfg["stage"] = stage
    selected = _selected_names(cfg, names)
    written: List[Path] = []
    for name in selected:
        spec = lookup_scorer(name)
        pattern = tool_inputs(spec, name)
        hits = expand_scoring_inputs(run_dir, pattern)
        if not hits:
            print(
                f"[scoring] skip {name}: no files matched {pattern!r} under {run_dir}",
                file=sys.stderr,
            )
            continue
        print(f"[scoring] {name} <- {', '.join(str(p) for p in hits)}")
        written.append(
            run_scorer(name, run_dir=run_dir, config=cfg, inputs=hits)
        )
    return written


def _selected_names(config: Dict[str, Any], names: Optional[Sequence[str]]) -> List[str]:
    if names is not None:
        wanted = [str(n).strip() for n in names if str(n).strip()]
    else:
        raw = config.get("scoring_tools")
        if raw in (None, "", False):
            wanted = iter_custom_scoring_names()
        elif isinstance(raw, (list, tuple)):
            wanted = [str(n).strip() for n in raw if str(n).strip()]
        else:
            wanted = [p.strip() for p in str(raw).replace(",", " ").split() if p.strip()]
    skip = {s.lower() for s in ("none", "off", "false", "0")}
    if any(w.lower() in skip for w in wanted):
        return []
    known = {n.lower(): n for n in iter_custom_scoring_names()}
    out: List[str] = []
    for item in wanted:
        matched = known.get(item.lower())
        if matched:
            out.append(matched)
        else:
            lookup_scorer(item)
            out.append(item)
    return out


def _try_python_callable(path: Path, name: str) -> Optional[Callable]:
    if path.suffix.lower() != ".py" or not path.is_file():
        return None
    try:
        spec = importlib.util.spec_from_file_location(
            f"samovar_custom_scorer_{name}", path
        )
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except Exception:
        return None
    for attr in ("score", "visualize", "visualise"):
        fn = getattr(module, attr, None)
        if callable(fn):
            return fn
    cls = getattr(module, "Scorer", None) or getattr(module, "Visualizer", None)
    if cls is None:
        return None
    try:
        inst = cls()
    except Exception:
        return None
    run = (
        getattr(inst, "score", None)
        or getattr(inst, "visualize", None)
        or getattr(inst, "run", None)
    )
    return run if callable(run) else None


def _run_cli_scorer(
    path: Path,
    inputs: Sequence[PathLike],
    output_dir: Path,
    config: Dict[str, Any],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [str(path)]
    for item in inputs:
        cmd.extend(["-i", str(item)])
    cmd.extend(["-o", str(output_dir)])
    extra = list(config.get("extra_argv") or extra_flags_argv(config.get("extra_flags")))
    cmd.extend(extra)
    if path.suffix.lower() == ".py":
        cmd = [sys.executable, *cmd]
    subprocess.run(cmd, check=True)


def load_scoring_config(path: Optional[PathLike]) -> Dict[str, Any]:
    if not path:
        return {}
    file = Path(path)
    if not file.is_file():
        return {}
    data = yaml.safe_load(file.read_text(encoding="utf-8")) or {}
    return data if isinstance(data, dict) else {}


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="samovar.scorers")
    sub = parser.add_subparsers(dest="command")
    run = sub.add_parser("run", help="Run imported scoring/viz tools on a run directory")
    add_output_dir_argument(run, required=True)
    run.add_argument("--config", default="", help="config_scoring.yaml from prepare")
    run.add_argument("--stage", default="", help="Pipeline checkpoint name (viz_initial, …)")
    run.add_argument(
        "--name",
        action="append",
        dest="names",
        default=None,
        help="Only these imported scoring tools (repeatable). Default: all custom.",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command != "run":
        build_parser().print_help()
        return 2
    cfg = load_scoring_config(args.config)
    cfg["stage"] = args.stage or cfg.get("stage") or ""
    run_dir = args.output_dir or cfg.get("output_dir")
    if not run_dir:
        print("Error: --output_dir is required", file=sys.stderr)
        return 1
    try:
        run_custom_scorers(
            run_dir,
            config=cfg,
            names=args.names,
            stage=str(cfg.get("stage") or ""),
        )
    except (MissingScorerError, subprocess.CalledProcessError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
