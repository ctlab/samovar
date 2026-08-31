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
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Optional, Sequence, Union

PathLike = Union[str, os.PathLike]

CHECKPOINT_STEPS = (
    "setup_reads",
    "qc_initial",
    "annotate_initial",
    "combine_initial",
    "viz_initial",
    "abundance_tables",
    "regenerate_tables",
    "score_regenerated_tables",
    "seed_genomes",
    "regenerate_reads",
    "sort_reads",
    "qc_generated",
    "annotate_regenerated",
    "combine_regenerated",
    "viz_regenerated",
    "reprofile",
    "viz_reprofiled",
)

# User-facing aliases for ``samovar prepare --startpoint/--endpoint``.
STEP_ALIASES = {
    "start": "setup_reads",
    "setup": "setup_reads",
    "reads": "setup_reads",
    "initial": "setup_reads",
    "qc": "qc_initial",
    "qc_initial": "qc_initial",
    "trim": "qc_initial",
    "annotate": "annotate_initial",
    "annotation": "annotate_initial",
    "combine": "combine_initial",
    "tables": "combine_initial",
    "table": "combine_initial",
    "generate_tables": "combine_initial",
    "annotation_tables": "combine_initial",
    "viz": "viz_initial",
    "plots": "viz_initial",
    "abundance": "abundance_tables",
    "abundance_tables": "abundance_tables",
    "otu": "abundance_tables",
    "otu_tables": "abundance_tables",
    "observed_abundance": "abundance_tables",
    "regenerate_tables": "regenerate_tables",
    "regen_tables": "regenerate_tables",
    "regenerated_tables": "regenerate_tables",
    "table_regen": "regenerate_tables",
    "score_tables": "score_regenerated_tables",
    "table_score": "score_regenerated_tables",
    "table_scoring": "score_regenerated_tables",
    "regenerated_tables_scoring": "score_regenerated_tables",
    "score_regenerated_tables": "score_regenerated_tables",
    "genomes": "seed_genomes",
    "seed": "seed_genomes",
    "regenerate": "regenerate_reads",
    "iss": "regenerate_reads",
    "simulation": "regenerate_reads",
    "sort": "sort_reads",
    "qc_generated": "qc_generated",
    "qc_regen": "qc_generated",
    "reannotate": "annotate_regenerated",
    "tables_regen": "combine_regenerated",
    "combine_regen": "combine_regenerated",
    "ml": "reprofile",
    "ensemble": "reprofile",
    "end": "viz_reprofiled",
}

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
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{name}.done.", suffix=".tmp", dir=str(dest)
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            handle.write(stamp + "\n")
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, marker)
    except Exception:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        raise
    if name in {"seed_genomes", "regenerate_reads"}:
        try:
            from samovar.genome_index import harvest_run_genomes

            harvest_run_genomes(output_dir)
        except Exception:
            pass
    try:
        from samovar.stage_report import write_stage_report

        write_stage_report(output_dir, name)
    except Exception:
        pass
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


def canonicalize_step(name: str) -> str:
    raw = str(name).strip()
    if not raw:
        raise ValueError("empty pipeline step name")
    key = raw.lower().replace("-", "_").replace(" ", "_")
    if key in CHECKPOINT_STEPS:
        return key
    alias = STEP_ALIASES.get(key)
    if alias:
        return alias
    known = ", ".join(CHECKPOINT_STEPS)
    extra = ", ".join(sorted(k for k in STEP_ALIASES if k not in CHECKPOINT_STEPS))
    raise ValueError(
        f"Unknown pipeline step {raw!r}. Canonical names: {known}. Aliases: {extra}."
    )


def resolve_window(
    start: Optional[str] = None,
    end: Optional[str] = None,
) -> tuple:
    start_raw = start if start not in (None, "") else os.environ.get("SAMOVAR_START")
    end_raw = end if end not in (None, "") else os.environ.get("SAMOVAR_END")
    start_s = canonicalize_step(start_raw or CHECKPOINT_STEPS[0])
    end_s = canonicalize_step(end_raw or CHECKPOINT_STEPS[-1])
    if CHECKPOINT_STEPS.index(start_s) > CHECKPOINT_STEPS.index(end_s):
        raise ValueError(f"startpoint {start_s} is after endpoint {end_s}")
    return start_s, end_s


def step_in_window(
    name: str,
    start: Optional[str] = None,
    end: Optional[str] = None,
) -> bool:
    start_s, end_s = resolve_window(start, end)
    idx = CHECKPOINT_STEPS.index(canonicalize_step(name))
    return CHECKPOINT_STEPS.index(start_s) <= idx <= CHECKPOINT_STEPS.index(end_s)


def needs_regen_early_exit(
    start: Optional[str] = None,
    end: Optional[str] = None,
) -> bool:
    """True when this run simulated reads and still has later regen stages."""
    if not step_in_window("regenerate_reads", start, end):
        return False
    later = CHECKPOINT_STEPS[CHECKPOINT_STEPS.index("regenerate_reads") + 1 :]
    return any(step_in_window(step, start, end) for step in later)


def _has_out_reports(directory: Path) -> bool:
    if not directory.is_dir():
        return False
    return any(path.is_file() for path in directory.rglob("*.out"))


def _has_csv_tables(directory: Path) -> bool:
    if not directory.is_dir():
        return False
    return any(path.is_file() and path.suffix == ".csv" for path in directory.glob("*.csv"))


def startpoint_gaps(
    output_dir: PathLike,
    start: Optional[str] = None,
    input_dir: Optional[str] = None,
) -> List[str]:
    """Describe missing on-disk inputs for ``start`` (empty list if ready)."""
    start_s, _ = resolve_window(start, CHECKPOINT_STEPS[-1])
    root = as_path(output_dir)
    from samovar.seqio import has_r1_reads

    src = input_dir if input_dir not in (None, "") else os.environ.get("SAMOVAR_INPUT_DIR")
    gaps: List[str] = []
    if start_s == "setup_reads":
        generate_sh = root / ".generate" / "generate.sh"
        if not (
            (src and has_r1_reads(src))
            or has_r1_reads(root / "initial")
            or generate_sh.is_file()
        ):
            gaps.append(
                "FASTQ under outdir/initial or --input_dir, or .generate/generate.sh"
            )
        return gaps
    if start_s == "qc_initial":
        if not has_r1_reads(root / "initial"):
            gaps.append("FASTQ samples under outdir/initial")
        return gaps
    if start_s == "annotate_initial":
        if not (
            has_r1_reads(root / "initial_trimmed") or has_r1_reads(root / "initial")
        ):
            gaps.append("FASTQ samples under outdir/initial_trimmed (or outdir/initial)")
        return gaps
    if start_s == "combine_initial":
        if not _has_out_reports(root / "initial_reports"):
            gaps.append("annotator *.out files under outdir/initial_reports")
        return gaps
    if start_s == "viz_initial":
        if not _has_csv_tables(root / "initial_annotations"):
            gaps.append("annotation CSVs under outdir/initial_annotations")
        return gaps
    if start_s == "abundance_tables":
        from samovar.abundance import collect_observed_abundance

        if not collect_observed_abundance(root):
            gaps.append(
                "abundance/OTU CSVs in outdir/initial_abundance (or outdir/*.csv / "
                "initial_annotations that convert to taxid + N_<sample>)"
            )
        return gaps
    if start_s == "regenerate_tables":
        from samovar.abundance import has_abundance_tables, observed_abundance_dir

        if not has_abundance_tables(observed_abundance_dir(root)):
            from samovar.abundance import collect_observed_abundance

            if not collect_observed_abundance(root):
                gaps.append(
                    "observed abundance CSVs under outdir/initial_abundance "
                    "(run abundance_tables, or drop OTU tables there)"
                )
        return gaps
    if start_s == "score_regenerated_tables":
        from samovar.abundance import (
            collect_observed_abundance,
            has_abundance_tables,
            observed_abundance_dir,
            regenerated_abundance_dir,
        )
        from samovar.table_scorers import load_tables_by_mode_from_run

        if not has_abundance_tables(observed_abundance_dir(root)) and not collect_observed_abundance(
            root
        ):
            gaps.append("observed abundance under outdir/initial_abundance")
        modes = load_tables_by_mode_from_run(root)
        if not modes and not has_abundance_tables(regenerated_abundance_dir(root)):
            gaps.append(
                "regenerated abundance CSVs under "
                "outdir/regenerated/.regenerated_abundance"
            )
        return gaps
    if start_s == "seed_genomes":
        return gaps
    if start_s == "regenerate_reads":
        from samovar.abundance import has_abundance_tables, regenerated_abundance_dir

        if not has_abundance_tables(regenerated_abundance_dir(root)):
            gaps.append(
                "regenerated abundance CSVs under "
                "outdir/regenerated/.regenerated_abundance "
                "(run regenerate_tables, or start earlier)"
            )
        return gaps
    if start_s == "sort_reads":
        if not has_r1_reads(root / "regenerated"):
            gaps.append("regenerated FASTQ under outdir/regenerated")
        return gaps
    if start_s == "qc_generated":
        if not has_r1_reads(root / "regenerated"):
            gaps.append("regenerated FASTQ under outdir/regenerated")
        return gaps
    if start_s == "annotate_regenerated":
        if not (
            has_r1_reads(root / "regenerated_trimmed")
            or has_r1_reads(root / "regenerated")
        ):
            gaps.append(
                "regenerated FASTQ under outdir/regenerated_trimmed (or outdir/regenerated)"
            )
        return gaps
    if start_s == "combine_regenerated":
        if not _has_out_reports(root / "regenerated_reports"):
            gaps.append("annotator *.out files under outdir/regenerated_reports")
        return gaps
    if start_s == "viz_regenerated":
        if not _has_csv_tables(root / "regenerated_annotations"):
            gaps.append("annotation CSVs under outdir/regenerated_annotations")
        return gaps
    if start_s == "reprofile":
        if not _has_csv_tables(root / "initial_annotations"):
            gaps.append("annotation CSVs under outdir/initial_annotations")
        regen_combined = root / "regenerated_annotations" / "combined_annotation_table.csv"
        if not regen_combined.is_file() and not _has_csv_tables(root / "regenerated_annotations"):
            gaps.append(
                "regenerated combined table "
                "(outdir/regenerated_annotations/combined_annotation_table.csv)"
            )
        return gaps
    if start_s == "viz_reprofiled":
        if not _has_csv_tables(root / "reprofiled_annotations"):
            gaps.append("annotation CSVs under outdir/reprofiled_annotations")
        return gaps
    return gaps


def check_startpoint(
    output_dir: PathLike,
    start: Optional[str] = None,
    end: Optional[str] = None,
    input_dir: Optional[str] = None,
) -> List[str]:
    """Validate the exec window and startpoint inputs. Returns error strings."""
    try:
        start_s, _end_s = resolve_window(start, end)
    except ValueError as exc:
        return [str(exc)]
    return startpoint_gaps(output_dir, start_s, input_dir=input_dir)


def run_checkup(
    output_dir: PathLike,
    start: Optional[str] = None,
    end: Optional[str] = None,
    input_dir: Optional[str] = None,
) -> int:
    """CLI-facing startpoint check. Exit 0 if ready, 1 if gaps, 2 if bad names."""
    try:
        start_s, end_s = resolve_window(start, end)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    gaps = startpoint_gaps(output_dir, start_s, input_dir=input_dir)
    root = as_path(output_dir)
    print(f"outdir: {root}")
    print(f"window: {start_s} .. {end_s}")
    if gaps:
        print(f"MISSING inputs for --start-point {start_s}:")
        for gap in gaps:
            print(f"  - {gap}")
        return 1
    print(f"OK: enough inputs to start at {start_s}")
    return 0


def _is_tmp_dir(path: Path) -> bool:
    return path.name in TMP_DIR_NAMES


def _is_tmp_file(name: str) -> bool:
    return "iss.tmp" in name


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
    pipe = sub.add_parser("pipeline", help="Print the exec script path (samovar.sh or samovar_vN.sh)")
    pipe.add_argument("output_dir")
    active = sub.add_parser(
        "active", help="Exit 0 if the named step is inside SAMOVAR_START/END"
    )
    active.add_argument("name")
    active.add_argument("--start", default=None)
    active.add_argument("--end", default=None)
    require = sub.add_parser(
        "require",
        help="Exit 1 if outdir is missing inputs for the startpoint",
    )
    require.add_argument("output_dir")
    require.add_argument("--start", default=None)
    require.add_argument("--end", default=None)
    require.add_argument("--input-dir", default=None)
    check = sub.add_parser(
        "checkup",
        help="Print whether outdir has the inputs required for a startpoint",
    )
    check.add_argument("output_dir", nargs="?", default=None)
    from samovar.paths import add_output_dir_argument

    add_output_dir_argument(check, dest="outdir_flag", default=None, required=False)
    check.add_argument("--start", "--startpoint", "--start-point", dest="start", default=None)
    check.add_argument("--end", "--endpoint", "--end-point", dest="end", default=None)
    check.add_argument("--input-dir", "--input_dir", dest="input_dir", default=None)
    regen = sub.add_parser(
        "needs-regen",
        help="Exit 0 if missing regenerated FASTQ should stop later stages",
    )
    regen.add_argument("--start", default=None)
    regen.add_argument("--end", default=None)
    steps = sub.add_parser("steps", help="Print canonical checkpoint names")
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
    if args.command == "pipeline":
        from samovar.add_annotator import active_pipeline_script

        print(active_pipeline_script(args.output_dir))
        return 0
    if args.command == "active":
        try:
            return 0 if step_in_window(args.name, args.start, args.end) else 1
        except ValueError as exc:
            print(exc, file=sys.stderr)
            return 2
    if args.command == "require":
        gaps = check_startpoint(
            args.output_dir,
            start=args.start,
            end=args.end,
            input_dir=args.input_dir,
        )
        if not gaps:
            return 0
        start_s = args.start or os.environ.get("SAMOVAR_START") or CHECKPOINT_STEPS[0]
        print(
            f"Error: startpoint {start_s!r} is missing required inputs in {args.output_dir}:",
            file=sys.stderr,
        )
        for gap in gaps:
            print(f"  - {gap}", file=sys.stderr)
        return 1
    if args.command == "checkup":
        return run_checkup(
            args.output_dir or args.outdir_flag or ".",
            start=args.start,
            end=args.end,
            input_dir=args.input_dir,
        )
    if args.command == "needs-regen":
        try:
            return 0 if needs_regen_early_exit(args.start, args.end) else 1
        except ValueError as exc:
            print(exc, file=sys.stderr)
            return 2
    if args.command == "steps":
        for name in CHECKPOINT_STEPS:
            print(name)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
