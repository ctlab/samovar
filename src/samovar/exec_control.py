"""Checkpoint markers and optional temporary-directory cleanup for ``samovar exec``.

Completed steps write ``$output_dir/.log/checkpoints/<name>.done``. A later
``samovar exec`` skips those steps unless ``--redo`` (or ``SAMOVAR_REDO=1``).
``--cleanup-tmp`` removes scratch dirs such as ``.tmp`` and ``.iss_full``.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence, Union

PathLike = Union[str, os.PathLike]

CHECKPOINT_STEPS = (
    "setup_reads",
    "annotate_initial",
    "combine_initial",
    "viz_initial",
    "seed_genomes",
    "regenerate_reads",
    "sort_reads",
    "annotate_regenerated",
    "combine_regenerated",
    "viz_regenerated",
    "reprofile",
    "viz_reprofiled",
)

CHECKPOINT_SUBDIR = Path(".log") / "checkpoints"
PROTECTED_TOPLEVEL = {".log", ".generate"}
TMP_DIR_NAMES = {".tmp", ".iss_full", ".combine_tmp"}


def as_path(path: PathLike) -> Path:
    return Path(path)


def checkpoint_dir(output_dir: PathLike) -> Path:
    return as_path(output_dir) / CHECKPOINT_SUBDIR


def redo_requested(explicit: Optional[bool] = None) -> bool:
    if explicit is not None:
        return bool(explicit)
    return os.environ.get("SAMOVAR_REDO", "0") not in {"", "0", "false", "False", "no"}


def should_skip(output_dir: PathLike, name: str, redo: Optional[bool] = None) -> bool:
    if redo_requested(redo):
        return False
    return (checkpoint_dir(output_dir) / f"{name}.done").is_file()


def mark_done(output_dir: PathLike, name: str) -> Path:
    dest = checkpoint_dir(output_dir)
    dest.mkdir(parents=True, exist_ok=True)
    marker = dest / f"{name}.done"
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    marker.write_text(stamp + "\n", encoding="utf-8")
    return marker


def clear_checkpoints(output_dir: PathLike, name: Optional[str] = None) -> List[Path]:
    dest = checkpoint_dir(output_dir)
    removed: List[Path] = []
    if name:
        marker = dest / f"{name}.done"
        if marker.is_file():
            marker.unlink()
            removed.append(marker)
        return removed
    if not dest.is_dir():
        return removed
    for marker in dest.glob("*.done"):
        marker.unlink()
        removed.append(marker)
    return removed


def listed_done(output_dir: PathLike) -> List[str]:
    dest = checkpoint_dir(output_dir)
    if not dest.is_dir():
        return []
    names = sorted(p.stem for p in dest.glob("*.done"))
    return names


def _is_tmp_dir(path: Path) -> bool:
    name = path.name
    if name in TMP_DIR_NAMES:
        return True
    if name.startswith("tmp_"):
        return True
    if name.endswith(".tmp"):
        return True
    return False


def _is_tmp_file(name: str) -> bool:
    return "iss.tmp" in name or name.endswith(".tmp")


def cleanup_tmp(output_dir: PathLike) -> List[str]:
    """Remove scratch directories and ISS temp files under ``output_dir``.

    Never deletes ``.log`` or ``.generate``. Paths outside the run root are ignored.
    """
    root = as_path(output_dir).resolve()
    removed: List[str] = []
    if not root.is_dir():
        return removed

    candidates: List[Path] = []
    for dirpath, dirnames, filenames in os.walk(root, topdown=True):
        current = Path(dirpath)
        rel = current.relative_to(root)
        if rel.parts and rel.parts[0] in PROTECTED_TOPLEVEL:
            dirnames[:] = []
            continue
        dirnames[:] = [d for d in dirnames if d not in PROTECTED_TOPLEVEL]
        for dirname in list(dirnames):
            child = current / dirname
            if _is_tmp_dir(child):
                candidates.append(child)
                dirnames.remove(dirname)
        for filename in filenames:
            if _is_tmp_file(filename):
                candidates.append(current / filename)

    for name in TMP_DIR_NAMES:
        hit = root / name
        if hit.exists():
            candidates.append(hit)

    seen = set()
    for path in sorted(candidates, key=lambda p: len(p.parts), reverse=True):
        try:
            resolved = path.resolve()
        except OSError:
            continue
        if resolved == root or resolved in seen:
            continue
        try:
            resolved.relative_to(root)
        except ValueError:
            continue
        seen.add(resolved)
        if resolved.is_dir():
            shutil.rmtree(resolved, ignore_errors=True)
            removed.append(str(resolved))
        elif resolved.is_file():
            try:
                resolved.unlink()
                removed.append(str(resolved))
            except OSError:
                pass
    return removed


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m samovar.exec_control",
        description="SamovaR exec checkpoints and temporary-directory cleanup",
    )
    sub = parser.add_subparsers(dest="command", required=True)
    skip = sub.add_parser("skip", help="Exit 0 if the named step should be skipped")
    skip.add_argument("output_dir")
    skip.add_argument("name")
    mark = sub.add_parser("mark", help="Record a completed checkpoint")
    mark.add_argument("output_dir")
    mark.add_argument("name")
    clear = sub.add_parser("clear", help="Remove checkpoint markers")
    clear.add_argument("output_dir")
    clear.add_argument("name", nargs="?")
    cleanup = sub.add_parser("cleanup", help="Remove .tmp / ISS scratch under output_dir")
    cleanup.add_argument("output_dir")
    listing = sub.add_parser("list", help="Print completed checkpoint names")
    listing.add_argument("output_dir")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(list(argv) if argv is not None else None)
    if args.command == "skip":
        return 0 if should_skip(args.output_dir, args.name) else 1
    if args.command == "mark":
        mark_done(args.output_dir, args.name)
        return 0
    if args.command == "clear":
        clear_checkpoints(args.output_dir, args.name)
        return 0
    if args.command == "cleanup":
        removed = cleanup_tmp(args.output_dir)
        for path in removed:
            print(path)
        return 0
    if args.command == "list":
        for name in listed_done(args.output_dir):
            print(name)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
