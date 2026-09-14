"""Merge two or more completed SamovaR runs into a new directory.

Modes:

* ``initial`` — union of initial-stage artifacts (including combined
  annotation tables). Downstream stages are left unmarked so ``samovar exec``
  rebuilds abundance, regeneration, and ML.
* ``regenerated`` — union of initial **and** regenerated artifacts (both
  combined annotation tables). Checkpoints through ``combine_regenerated`` are
  marked so exec continues at visualization and reprofiling.
"""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple

import pandas as pd
import yaml

from samovar.abundance import (
    load_abundance_dir,
    n_sample_columns,
    normalize_abundance_table,
    write_abundance_dir,
)
from samovar.add_annotator import existing_run_names, write_active_pipeline
from samovar.exec_control import (
    CHECKPOINT_STEPS,
    as_path,
    clear_checkpoints,
    mark_done,
)
from samovar.paths import absolute_path, add_output_dir_argument
from samovar.seqio import list_fastq_samples

PathLike = str | Path

SAMPLE_NAME_ERROR = (
    "Sample names need to be distinct. Names solutions is not yet implemented"
)
COMBINED_PREFIX = "combined_annotation_table"
SKIP_DIR_NAMES = {
    ".snakemake",
    ".tmp",
    ".cache",
    ".iss_full",
    ".combine_tmp",
    "__pycache__",
}
SKIP_FILE_NAMES = {".process_annotations.done"}
SCAFFOLD_RELS = (".log", ".generate", ".database", ".hydra")
INITIAL_RELS = (
    "initial",
    "initial_trimmed",
    "initial_reports",
    "initial_annotations",
)
REGENERATED_RELS = (
    "regenerated",
    "regenerated_trimmed",
    "regenerated_reports",
    "regenerated_annotations",
    "genomes",
)
ABUNDANCE_RELS = {
    "initial_abundance",
    str(Path("regenerated") / ".regenerated_abundance"),
}
KEEP_EXISTING_TOP = {"genomes", ".database", ".generate", ".hydra"}
REWRITE_SUFFIXES = {
    ".yaml",
    ".yml",
    ".sh",
    ".env",
    ".json",
    ".txt",
    ".csv",
}
MODE_ALIASES = {
    "initial": "initial",
    "init": "initial",
    "regenerated": "regenerated",
    "regen": "regenerated",
    "regenerate": "regenerated",
}


class MergeError(ValueError):
    """Invalid merge inputs or colliding sample names."""


def canonicalize_mode(name: str) -> str:
    key = str(name or "").strip().lower().replace("-", "_")
    mapped = MODE_ALIASES.get(key)
    if mapped is None:
        known = ", ".join(sorted(set(MODE_ALIASES.values())))
        raise MergeError(f"Unknown merge mode {name!r}. Use {known}.")
    return mapped


def _is_combined_table(name: str) -> bool:
    return Path(name).name.startswith(COMBINED_PREFIX)


def _annotation_sample_names(ann_dir: Path) -> List[str]:
    if not ann_dir.is_dir():
        return []
    names: List[str] = []
    for path in sorted(ann_dir.glob("*.csv")):
        if _is_combined_table(path.name):
            continue
        if path.name.endswith(".annotation.csv"):
            names.append(path.name[: -len(".annotation.csv")])
        else:
            names.append(path.stem)
    return names


def _abundance_sample_names(folder: Path) -> List[str]:
    names: List[str] = []
    for table in load_abundance_dir(folder).values():
        for col in n_sample_columns(table):
            names.append(str(col)[2:] if str(col).startswith("N_") else str(col))
    return names


def _combined_sample_names(path: Path) -> List[str]:
    if not path.is_file():
        return []
    try:
        frame = pd.read_csv(path, usecols=lambda c: str(c).lower() == "sample")
    except (ValueError, pd.errors.EmptyDataError, OSError):
        return []
    if frame.empty or "sample" not in {str(c).lower() for c in frame.columns}:
        return []
    col = next(c for c in frame.columns if str(c).lower() == "sample")
    return [str(v) for v in frame[col].dropna().unique() if str(v).strip()]


def sample_names_for_run(output_dir: PathLike, mode: str) -> List[str]:
    root = as_path(output_dir)
    mode_s = canonicalize_mode(mode)
    names: List[str] = []
    names.extend(list_fastq_samples(root / "initial"))
    names.extend(list_fastq_samples(root / "initial_trimmed"))
    names.extend(_annotation_sample_names(root / "initial_annotations"))
    names.extend(_abundance_sample_names(root / "initial_abundance"))
    names.extend(
        _combined_sample_names(
            root / "initial_annotations" / "combined_annotation_table.csv"
        )
    )
    if mode_s == "regenerated":
        names.extend(list_fastq_samples(root / "regenerated"))
        names.extend(list_fastq_samples(root / "regenerated_trimmed"))
        names.extend(_annotation_sample_names(root / "regenerated_annotations"))
        names.extend(
            _abundance_sample_names(
                root / "regenerated" / ".regenerated_abundance"
            )
        )
        names.extend(
            _combined_sample_names(
                root / "regenerated_annotations" / "combined_annotation_table.csv"
            )
        )
    seen: Set[str] = set()
    out: List[str] = []
    for name in names:
        token = str(name).strip()
        if not token or token in seen:
            continue
        seen.add(token)
        out.append(token)
    return out


def assert_distinct_sample_names(sources: Sequence[Path], mode: str) -> List[str]:
    owner: Dict[str, Path] = {}
    ordered: List[str] = []
    for src in sources:
        for name in sample_names_for_run(src, mode):
            prev = owner.get(name)
            if prev is not None and prev != src:
                raise MergeError(SAMPLE_NAME_ERROR)
            if name not in owner:
                owner[name] = src
                ordered.append(name)
    return ordered


def assert_matching_annotators(sources: Sequence[Path]) -> List[str]:
    sets = [tuple(existing_run_names(src)) for src in sources]
    primary = list(sets[0])
    for names in sets[1:]:
        if set(names) != set(primary):
            raise MergeError(
                "cannot merge runs with different annotator sets "
                f"(have {sorted(set(primary))} vs {sorted(set(names))})"
            )
    return primary


def _rels_for_mode(mode: str) -> List[str]:
    rels = list(INITIAL_RELS)
    if mode == "regenerated":
        rels.extend(REGENERATED_RELS)
        rels.append("initial_abundance")
    return rels


def _keep_through(step: str) -> Set[str]:
    idx = CHECKPOINT_STEPS.index(step)
    return set(CHECKPOINT_STEPS[: idx + 1])


def _is_abundance_rel(rel: Path) -> bool:
    posix = rel.as_posix()
    if posix in ABUNDANCE_RELS:
        return True
    return rel.parent.as_posix() in ABUNDANCE_RELS


def _should_skip_dir(name: str) -> bool:
    return name in SKIP_DIR_NAMES or name in SKIP_FILE_NAMES


def _copy_file(src: Path, dest: Path, *, allow_existing: bool) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        if allow_existing:
            return
        raise MergeError(SAMPLE_NAME_ERROR)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_symlink():
        shutil.copy2(src, dest, follow_symlinks=False)
        return
    shutil.copy2(src, dest)


def _merge_abundance_frames(frames: Sequence[pd.DataFrame]) -> pd.DataFrame:
    out: Optional[pd.DataFrame] = None
    for raw in frames:
        frame = normalize_abundance_table(raw)
        if frame.empty:
            continue
        if out is None:
            out = frame
            continue
        overlap = set(n_sample_columns(out)) & set(n_sample_columns(frame))
        if overlap:
            raise MergeError(SAMPLE_NAME_ERROR)
        out = out.merge(frame, on="taxid", how="outer")
    if out is None:
        return pd.DataFrame(columns=["taxid"])
    for col in n_sample_columns(out):
        out[col] = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
    return out.reset_index(drop=True)


def _concat_tables(paths: Sequence[Path]) -> Optional[pd.DataFrame]:
    frames: List[pd.DataFrame] = []
    for path in paths:
        try:
            frame = pd.read_csv(path)
        except (pd.errors.EmptyDataError, OSError):
            continue
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return None
    return pd.concat(frames, ignore_index=True, sort=False)


def _rewrite_paths(root: Path, sources: Sequence[Path], dest: Path) -> None:
    replacements: List[Tuple[str, str]] = []
    dest_s = str(dest.resolve())
    for src in sorted((p.resolve() for p in sources), key=lambda p: len(str(p)), reverse=True):
        text = str(src)
        if text and text != dest_s:
            replacements.append((text, dest_s))
    if not replacements:
        return
    for folder in (root / ".log", root / ".hydra", root / ".generate"):
        if not folder.is_dir():
            continue
        for path in folder.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() not in REWRITE_SUFFIXES and path.name not in {
                "active_pipeline",
                "window.env",
                "samovar.sh",
            }:
                if not path.name.endswith(".sh"):
                    continue
            try:
                body = path.read_text(encoding="utf-8")
            except (OSError, UnicodeDecodeError):
                continue
            updated = body
            for old, new in replacements:
                updated = updated.replace(old, new)
            if updated != body:
                path.write_text(updated, encoding="utf-8")


def _copy_scaffold(primary: Path, dest: Path) -> None:
    for rel in SCAFFOLD_RELS:
        src = primary / rel
        if src.exists():
            shutil.copytree(src, dest / rel, dirs_exist_ok=True, symlinks=True)
    ckpt = dest / ".log" / "checkpoints"
    if ckpt.is_dir():
        shutil.rmtree(ckpt)


def _iter_files(folder: Path) -> Iterable[Path]:
    for path in folder.rglob("*"):
        if not path.is_file() and not path.is_symlink():
            continue
        rel_parts = path.relative_to(folder).parts
        if any(_should_skip_dir(part) for part in rel_parts):
            continue
        if path.name in SKIP_FILE_NAMES:
            continue
        yield path


def _merge_trees(
    sources: Sequence[Path],
    dest: Path,
    rels: Sequence[str],
) -> Tuple[Dict[Path, List[Path]], Dict[Path, List[pd.DataFrame]]]:
    combined: Dict[Path, List[Path]] = {}
    abundance: Dict[Path, List[pd.DataFrame]] = {}
    for rel in rels:
        for src in sources:
            src_root = src / rel
            if not src_root.exists():
                continue
            if src_root.is_file():
                _copy_file(
                    src_root,
                    dest / rel,
                    allow_existing=rel.split("/")[0] in KEEP_EXISTING_TOP,
                )
                continue
            for path in _iter_files(src_root):
                rel_path = Path(rel) / path.relative_to(src_root)
                if _is_combined_table(path.name):
                    combined.setdefault(dest / rel_path, []).append(path)
                    continue
                if _is_abundance_rel(rel_path):
                    if path.suffix == ".csv":
                        try:
                            frame = pd.read_csv(path)
                        except (pd.errors.EmptyDataError, OSError):
                            continue
                        abundance.setdefault(dest / rel_path, []).append(frame)
                    else:
                        _copy_file(
                            path, dest / rel_path, allow_existing=True
                        )
                    continue
                allow = rel_path.parts[0] in KEEP_EXISTING_TOP
                _copy_file(path, dest / rel_path, allow_existing=allow)
    return combined, abundance


def _write_combined(groups: Dict[Path, List[Path]]) -> None:
    for dest, paths in groups.items():
        table = _concat_tables(paths)
        if table is None:
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        table.to_csv(dest, index=False)


def _write_abundance(groups: Dict[Path, List[pd.DataFrame]]) -> None:
    by_dir: Dict[Path, Dict[str, List[pd.DataFrame]]] = {}
    for dest, frames in groups.items():
        by_dir.setdefault(dest.parent, {}).setdefault(dest.stem, []).extend(frames)
    for folder, tables in by_dir.items():
        merged = {
            stem: _merge_abundance_frames(frames) for stem, frames in tables.items()
        }
        write_abundance_dir(folder, merged)


def _write_checkpoints(dest: Path, mode: str) -> List[str]:
    keep = _keep_through(
        "combine_initial" if mode == "initial" else "combine_regenerated"
    )
    clear_checkpoints(dest)
    marked: List[str] = []
    for name in CHECKPOINT_STEPS:
        if name in keep:
            mark_done(dest, name)
            marked.append(name)
    return marked


def merge_runs(
    sources: Sequence[PathLike],
    output_dir: PathLike,
    mode: str,
    *,
    force: bool = False,
) -> Dict[str, Any]:
    mode_s = canonicalize_mode(mode)
    runs = [as_path(s).expanduser().resolve() for s in sources]
    if len(runs) < 2:
        raise MergeError("samovar merge needs two or more run directories")
    dest = as_path(output_dir).expanduser().resolve()
    for src in runs:
        if not src.is_dir():
            raise MergeError(f"run directory not found: {src}")
        cfg = src / ".log" / "configs" / "config_init.yaml"
        if not cfg.is_file():
            raise MergeError(
                f"{src} has no .log/configs/config_init.yaml "
                "(prepare/exec a source run first)"
            )
        if src == dest:
            raise MergeError("merge --output_dir must not be one of the source runs")
    samples = assert_distinct_sample_names(runs, mode_s)
    annotators = assert_matching_annotators(runs)
    if dest.exists() and any(dest.iterdir()):
        if not force:
            raise MergeError(
                f"destination {dest} is not empty (pass --force to replace it)"
            )
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)

    primary = runs[0]
    _copy_scaffold(primary, dest)
    combined, abundance = _merge_trees(runs, dest, _rels_for_mode(mode_s))
    _write_combined(combined)
    _write_abundance(abundance)
    _rewrite_paths(dest, runs, dest)
    pipe = dest / ".log" / "samovar.sh"
    if pipe.is_file():
        write_active_pipeline(dest, "samovar.sh")
    marked = _write_checkpoints(dest, mode_s)
    meta = {
        "mode": mode_s,
        "sources": [str(p) for p in runs],
        "output_dir": str(dest),
        "samples": samples,
        "annotators": annotators,
        "checkpoints": marked,
        "next_step": (
            "viz_initial" if mode_s == "initial" else "viz_regenerated"
        ),
    }
    dest_log = dest / ".log"
    dest_log.mkdir(parents=True, exist_ok=True)
    (dest_log / "merge.yaml").write_text(
        yaml.safe_dump(meta, sort_keys=False), encoding="utf-8"
    )
    return meta


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="samovar merge",
        description=(
            "Combine two or more SamovaR run directories. "
            "Sample names must be distinct. "
            "Mode 'initial' unions initial-stage annotations; "
            "mode 'regenerated' also unions regenerated annotations."
        ),
    )
    parser.add_argument(
        "--mode",
        "-m",
        required=True,
        help="initial | regenerated",
    )
    add_output_dir_argument(parser, required=True)
    parser.add_argument(
        "--force",
        action="store_true",
        help="Replace the destination if it already exists",
    )
    parser.add_argument(
        "runs",
        nargs="+",
        help="Source run directories (at least two)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(list(argv) if argv is not None else None)
    try:
        meta = merge_runs(
            args.runs,
            args.output_dir,
            args.mode,
            force=bool(args.force),
        )
    except MergeError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(
        f"merge ({meta['mode']}): {len(meta['sources'])} runs -> "
        f"{absolute_path(args.output_dir)}; "
        f"next {meta['next_step']}; samples {', '.join(meta['samples']) or '(none)'}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
