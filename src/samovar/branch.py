"""Copy a run up to a start-point and continue prepare from that step.

``samovar branch --start-point regenerate_tables SRC DST`` copies artifacts of
completed steps before ``regenerate_tables``. A later ``samovar prepare`` in
``DST`` keeps those files, merges Hydra prepare args, rewrites ``samovar.sh``
with that window, and refuses to change the initial annotator set.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from argparse import Namespace
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import yaml

from samovar.add_annotator import (
    cli_annotator_keys,
    existing_run_names,
    load_prepare_hydra_args,
    log_dir,
    write_active_pipeline,
)
from samovar.exec_control import (
    CHECKPOINT_STEPS,
    as_path,
    canonicalize_step,
    check_startpoint,
    clear_checkpoints,
)
from samovar.paths import absolute_path, add_output_dir_argument

BRANCH_META = Path(".log") / "branch.yaml"
SKIP_TOPLEVEL = {
    ".snakemake",
    ".cache",
    ".tmp",
    ".iss_full",
    "multiqc_samovar",
    "multiqc",
}
UNSPECIFIED_FALSE = {
    "use_test_genomes",
    "skip_export",
    "to_mpa",
    "to_cami",
    "to_kraken2",
    "to_abundance",
    "add_annotator",
}


def branch_meta_path(output_dir) -> Path:
    return as_path(output_dir) / BRANCH_META


def is_branch_dir(output_dir) -> bool:
    return branch_meta_path(output_dir).is_file()


def load_branch_meta(output_dir) -> Dict[str, Any]:
    path = branch_meta_path(output_dir)
    if not path.is_file():
        return {}
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    return data if isinstance(data, dict) else {}


def write_branch_meta(output_dir, payload: Dict[str, Any]) -> Path:
    dest = branch_meta_path(output_dir)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(yaml.safe_dump(payload, sort_keys=False), encoding="utf-8")
    return dest


def annotator_names(output_dir) -> List[str]:
    return list(existing_run_names(output_dir))


def _run_names_from_cli(args: Any) -> List[str]:
    names: List[str] = []
    for key in cli_annotator_keys(args):
        names.append(key[4:] if str(key).startswith("cmd_") else str(key))
    seen = set()
    out: List[str] = []
    for name in names:
        if name not in seen:
            seen.add(name)
            out.append(name)
    return out


def reject_annotator_change(output_dir, args: Any) -> None:
    """Branches cannot swap the initial classifier set (use add-annotator on the source)."""
    frozen = load_branch_meta(output_dir).get("annotators") or annotator_names(output_dir)
    frozen_set = set(str(n) for n in frozen if n)
    offered = _run_names_from_cli(args)
    if not offered:
        return
    if not frozen_set:
        return
    if set(offered) != frozen_set:
        raise ValueError(
            "cannot change initial annotators on a branched run "
            f"(have {sorted(frozen_set)}, got {sorted(offered)}). "
            "Use samovar prepare --add-annotator on the original directory."
        )


def paths_produced_before(start: str) -> List[str]:
    """Relative output dirs written by checkpoints strictly before ``start``."""
    from samovar.stage_report import STAGE_DIRS

    start_s = canonicalize_step(start)
    start_i = CHECKPOINT_STEPS.index(start_s)
    seen: Set[str] = set()
    ordered: List[str] = []
    for step in CHECKPOINT_STEPS[:start_i]:
        for rel in STAGE_DIRS.get(step, ()):
            if rel not in seen:
                seen.add(rel)
                ordered.append(rel)
    # If the parent tree is copied, drop redundant children.
    collapsed: List[str] = []
    for rel in sorted(ordered, key=lambda p: (p.count("/"), p)):
        if any(
            rel != parent and (rel.startswith(parent + "/") or rel.startswith(parent + "\\"))
            for parent in collapsed
        ):
            continue
        collapsed.append(rel)
    return collapsed


def _ignore_scratch(directory: str, names: List[str]) -> Set[str]:
    skip = set()
    for name in names:
        if name in SKIP_TOPLEVEL or name.endswith(".iss_full"):
            skip.add(name)
        if name == ".process_annotations.done":
            skip.add(name)
    return skip


def _copy_tree(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_file():
        shutil.copy2(src, dest)
        return
    if dest.exists():
        shutil.copytree(
            src,
            dest,
            dirs_exist_ok=True,
            ignore=_ignore_scratch,
            symlinks=True,
        )
        return
    shutil.copytree(src, dest, ignore=_ignore_scratch, symlinks=True)


def copy_run_until(
    source,
    dest,
    start: str,
    *,
    force: bool = False,
) -> Dict[str, Any]:
    start_s = canonicalize_step(start)
    src = as_path(source).resolve()
    dst = as_path(dest).resolve()
    if not src.is_dir():
        raise ValueError(f"branch source is not a directory: {src}")
    if src == dst:
        raise ValueError("branch source and destination must differ")
    gaps = check_startpoint(src, start=start_s)
    if gaps:
        raise ValueError(
            f"source is missing inputs for --start-point {start_s}: " + "; ".join(gaps)
        )
    if dst.exists() and any(dst.iterdir()):
        if not force:
            raise ValueError(
                f"destination {dst} is not empty (pass --force to replace it)"
            )
        shutil.rmtree(dst)
    dst.mkdir(parents=True, exist_ok=True)

    copied: List[str] = []
    always = (".hydra", ".generate", ".database", "reads")
    for rel in always:
        piece = src / rel
        if piece.exists():
            _copy_tree(piece, dst / rel)
            copied.append(rel)

    log_src = src / ".log"
    if log_src.is_dir():
        _copy_tree(log_src, dst / ".log")
        ckpt = dst / ".log" / "checkpoints"
        if ckpt.is_dir():
            start_i = CHECKPOINT_STEPS.index(start_s)
            keep = set(CHECKPOINT_STEPS[:start_i])
            for marker in list(ckpt.glob("*.done")):
                if marker.stem not in keep:
                    marker.unlink()
        copied.append(".log")

    for rel in paths_produced_before(start_s):
        piece = src / rel
        if piece.exists():
            _copy_tree(piece, dst / rel)
            copied.append(rel)

    meta = {
        "source": str(src),
        "startpoint": start_s,
        "annotators": annotator_names(src),
        "copied": copied,
    }
    write_branch_meta(dst, meta)
    return meta


def _overlay_cli(merged: Dict[str, Any], args: Any) -> Dict[str, Any]:
    data = vars(args) if not isinstance(args, dict) else dict(args)
    for key, value in data.items():
        if key in ("add_annotator",):
            continue
        if str(key).startswith("cmd_") or key in {"kraken2", "kaiju", "dummy"}:
            continue
        if value is None:
            continue
        if key in UNSPECIFIED_FALSE and value is False:
            continue
        merged[key] = value
    return merged


def maybe_merge_branch_prepare(args: Any) -> Tuple[Any, bool]:
    """Merge Hydra + CLI when dest is a branch (or startpoint on an existing run)."""
    out_dir = absolute_path(getattr(args, "output_dir", None) or ".")
    branched = is_branch_dir(out_dir)
    start_cli = getattr(args, "startpoint", None)
    hydra = load_prepare_hydra_args(out_dir)
    if getattr(args, "add_annotator", False) and branched:
        raise ValueError(
            "cannot --add-annotator on a branched directory; "
            "add classifiers on the original run"
        )
    if not branched and not (hydra and start_cli):
        return args, False
    if not hydra and not (log_dir(out_dir) / "samovar.sh").is_file():
        if branched:
            raise ValueError(f"branched run {out_dir} has no Hydra prepare snapshot")
        return args, False
    if not branched:
        # Existing startpoint prepare on a full run: rewrite the window, keep files.
        pass
    reject_annotator_change(out_dir, args)
    merged = dict(hydra)
    merged = _overlay_cli(merged, args)
    meta = load_branch_meta(out_dir)
    start_s = start_cli or meta.get("startpoint") or merged.get("startpoint")
    if start_s:
        merged["startpoint"] = canonicalize_step(str(start_s))
    merged["output_dir"] = out_dir
    merged.pop("add_annotator", None)
    ns = Namespace(**merged)
    ns.add_annotator = False
    ns.output_dir = out_dir
    if not getattr(ns, "input_dir", None) and not getattr(ns, "input_config", None):
        initial = as_path(out_dir) / "initial"
        if initial.is_dir():
            ns.input_dir = str(initial)
    return ns, branched


def finish_branch_prepare(output_dir, args: Any) -> None:
    start_s = canonicalize_step(
        getattr(args, "startpoint", None) or load_branch_meta(output_dir).get("startpoint")
        or CHECKPOINT_STEPS[0]
    )
    write_active_pipeline(output_dir, "samovar.sh")
    start_i = CHECKPOINT_STEPS.index(start_s)
    for name in CHECKPOINT_STEPS[start_i:]:
        clear_checkpoints(output_dir, name)
    meta = load_branch_meta(output_dir) or {}
    if meta:
        meta["startpoint"] = start_s
        meta["annotators"] = annotator_names(output_dir) or meta.get("annotators") or []
        write_branch_meta(output_dir, meta)


def _resolve_paths(args: argparse.Namespace) -> Tuple[Path, Path]:
    positionals = [Path(p) for p in (args.paths or [])]
    source = args.source
    dest = args.output_dir
    if len(positionals) == 2:
        source = source or str(positionals[0])
        dest = dest or str(positionals[1])
    elif len(positionals) == 1:
        if dest:
            source = source or str(positionals[0])
        else:
            dest = str(positionals[0])
    if not source:
        source = "."
    if not dest:
        raise ValueError("destination is required (--output_dir / --directory, or SRC DST)")
    return Path(source), Path(dest)


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        prog="samovar branch",
        description=(
            "Copy a SamovaR run up to --start-point. Then samovar prepare in the "
            "copy keeps those artifacts and rewrites the pipeline from that step."
        ),
    )
    parser.add_argument(
        "--start-point",
        "--startpoint",
        dest="startpoint",
        required=True,
        help="First step the copy will re-run (same names as prepare --startpoint)",
    )
    parser.add_argument(
        "--from",
        "--source",
        dest="source",
        default=None,
        help="Source run directory (default: first positional or .)",
    )
    add_output_dir_argument(
        parser,
        required=False,
        default=None,
        help="Destination run directory (or pass SRC DST)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace the destination if it already exists",
    )
    parser.add_argument("paths", nargs="*", help="Optional SRC DST")
    args = parser.parse_args(list(argv) if argv is not None else None)
    try:
        source, dest = _resolve_paths(args)
        meta = copy_run_until(source, dest, args.startpoint, force=args.force)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(
        f"branch: {meta['source']} -> {absolute_path(dest)}; "
        f"start-point {meta['startpoint']}; "
        f"annotators {', '.join(meta.get('annotators') or []) or '(none)'}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
