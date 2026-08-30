"""Read QC: metagenome FASTQ → trimmed metagenome FASTQ.

Default is identity (reads are already trimmed). Imported tools use
``samovar tools import --type QC``. Hybrid runs may pick a different tool
per filename postfix (``illumina``, ``ont``, …).
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple

from samovar.main_config import (
    flags_target_matches,
    imported_flags_for_names,
    iter_tools,
    merge_flag_strings,
    parse_tool_entry,
    tool_path,
)
from samovar.paths import load_config
from samovar.reads_generators import FILENAME_EXTRA_IDS
from samovar.seqio import (
    as_path,
    find_fastq_mate,
    is_fastq_name,
    iter_fastq_records,
    open_text,
    r2_for_r1,
    sample_name_from_r1,
)
from samovar.table_regenerators import extra_flags_argv

QC_GROUP = "qc"
IDENTITY_NAMES = frozenset({"", "none", "off", "false", "0", "identity", "trim"})
QC_FLAG_GROUPS = ("qc", "quality", "trim", "trimming")


class MissingQCError(ValueError):
    """``tools.<name>`` is missing or is not ``--type QC``."""


def is_identity_qc(name: Optional[str]) -> bool:
    return str(name or "").strip().lower() in IDENTITY_NAMES


def normalize_qc_name(name: Optional[str]) -> str:
    token = str(name or "").strip()
    if is_identity_qc(token):
        return ""
    return token


def parse_qc_postfix(raw: Any) -> Dict[str, str]:
    """``illumina:gc_filter`` tokens or ``{illumina: name}``."""
    out: Dict[str, str] = {}
    if raw in (None, False, ""):
        return out
    if isinstance(raw, dict):
        for key, value in raw.items():
            postfix = str(key or "").strip().lower()
            tool = normalize_qc_name(value)
            if postfix:
                out[postfix] = tool
        return out
    items: List[Any]
    if isinstance(raw, (list, tuple)):
        items = list(raw)
    else:
        items = [raw]
    for item in items:
        text = str(item or "").strip()
        if not text:
            continue
        if ":" in text:
            postfix, tool = text.split(":", 1)
        elif "=" in text:
            postfix, tool = text.split("=", 1)
        else:
            continue
        postfix = postfix.strip().lower()
        if postfix:
            out[postfix] = normalize_qc_name(tool)
    return out


def flags_apply_to_qc(target: str, *names: Optional[str]) -> bool:
    return flags_target_matches(
        target,
        *[str(n) for n in names if n],
        groups=QC_FLAG_GROUPS,
    )


def lookup_qc(name: str) -> list:
    key = str(name or "").strip()
    if not key:
        raise MissingQCError(
            "Empty QC name. Import a tool with "
            "`samovar tools import -n NAME --type QC`."
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
        raise MissingQCError(
            f"QC {key!r} is not in the main install config. "
            "Register it with `samovar tools import -n "
            f"{key} --exec-path /path/to/script.py --type QC` "
            "before prepare."
        )
    parsed = parse_tool_entry(spec, matched)
    group = str(parsed[3] or "").strip()
    if group != QC_GROUP:
        raise MissingQCError(
            f"tools.{matched} has group {group!r}, expected {QC_GROUP!r}. "
            "Re-import with --type QC."
        )
    path = tool_path(parsed, matched)
    if not path:
        raise MissingQCError(
            f"tools.{matched} has an empty path. Re-import with --exec-path."
        )
    return parsed


def require_known_qc(name: Optional[str]) -> str:
    token = normalize_qc_name(name)
    if not token:
        return ""
    lookup_qc(token)
    return token


def attach_qc_flags(name: Optional[str], config: Dict[str, Any]) -> Dict[str, Any]:
    cfg = dict(config or {})
    token = normalize_qc_name(name or cfg.get("qc"))
    tools = iter_tools(load_config())
    imported = imported_flags_for_names(tools, token) if token else ""
    named = ""
    if token:
        named = str((cfg.get("qc_tool_flags") or {}).get(token) or "")
    cfg["extra_flags"] = merge_flag_strings(
        imported,
        cfg.get("qc_flags"),
        cfg.get("extra_flags"),
        named,
    )
    cfg["extra_argv"] = extra_flags_argv(cfg.get("extra_flags"))
    _apply_gc_flag_tokens(cfg)
    cfg["qc"] = token
    return cfg


def _apply_gc_flag_tokens(cfg: Dict[str, Any]) -> None:
    argv = list(cfg.get("extra_argv") or [])
    i = 0
    while i < len(argv):
        tok = str(argv[i])
        if tok in {"--min-gc", "--min_gc", "--gc-min"} and i + 1 < len(argv):
            cfg["min_gc"] = float(argv[i + 1])
            i += 2
            continue
        if tok in {"--max-gc", "--max_gc", "--gc-max"} and i + 1 < len(argv):
            cfg["max_gc"] = float(argv[i + 1])
            i += 2
            continue
        i += 1


def postfix_from_sample(sample: str) -> str:
    text = str(sample or "")
    parts = text.split("_")
    extras = {str(x).lower() for x in FILENAME_EXTRA_IDS}
    for token in reversed(parts):
        low = token.lower()
        if low in extras:
            return low
    return ""


def qc_name_for_sample(
    sample: str,
    *,
    stage_qc: Optional[str] = None,
    postfix_map: Optional[Dict[str, str]] = None,
) -> str:
    postfix = postfix_from_sample(sample)
    mapping = postfix_map or {}
    if postfix and postfix in mapping:
        return normalize_qc_name(mapping[postfix])
    return normalize_qc_name(stage_qc)


def list_all_r1_files(directory: Path) -> List[Path]:
    """Every R1 FASTQ, including hybrid ``_{tech}`` copies dropped for annotators."""
    folder = as_path(directory)
    if not folder.is_dir():
        return []
    found: List[Path] = []
    for path in sorted(folder.iterdir()):
        if not path.is_file() or not is_fastq_name(path.name):
            continue
        try:
            sample_name_from_r1(path)
        except ValueError:
            continue
        found.append(path)
    return found


def write_fastq_records(path: Path, records: Sequence[Tuple[str, str, str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        for header, seq, plus, qual in records:
            for line in (header, seq, plus, qual):
                handle.write(line if str(line).endswith("\n") else str(line) + "\n")


def count_fastq_records(path: Path) -> int:
    if not path.is_file():
        return 0
    return sum(1 for _ in iter_fastq_records(path))


def gc_fraction(seq: str) -> float:
    bases = [c for c in str(seq).upper() if c in "ACGT"]
    if not bases:
        return 0.5
    return (bases.count("G") + bases.count("C")) / len(bases)


def _try_python_trim(path: Path, name: str) -> Optional[Callable]:
    if path.suffix.lower() != ".py" or not path.is_file():
        return None
    try:
        spec = importlib.util.spec_from_file_location(f"samovar_custom_qc_{name}", path)
        if spec is None or spec.loader is None:
            return None
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
    except Exception:
        return None
    fn = getattr(module, "trim", None)
    if callable(fn):
        return fn
    cls = getattr(module, "QC", None) or getattr(module, "ReadQC", None)
    if cls is None:
        return None
    try:
        inst = cls()
    except Exception:
        return None
    run = getattr(inst, "trim", None) or getattr(inst, "run", None)
    return run if callable(run) else None


def _link_or_copy_file(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() or dest.is_symlink():
        dest.unlink()
    try:
        dest.symlink_to(src.resolve())
    except OSError:
        import shutil

        shutil.copy2(src, dest)


def identity_pair(r1: Path, r2: Optional[Path], dest_r1: Path, dest_r2: Optional[Path]) -> List[str]:
    _link_or_copy_file(r1, dest_r1)
    if r2 is not None and r2.exists() and dest_r2 is not None:
        _link_or_copy_file(r2, dest_r2)
        return [str(dest_r1), str(dest_r2)]
    return [str(dest_r1)]


def run_qc_pair(
    name: str,
    r1: Path,
    r2: Optional[Path],
    dest_r1: Path,
    dest_r2: Optional[Path],
    config: Dict[str, Any],
) -> List[str]:
    token = normalize_qc_name(name)
    if not token:
        return identity_pair(r1, r2, dest_r1, dest_r2)
    spec = lookup_qc(token)
    path = Path(tool_path(spec, token)).expanduser()
    cfg = attach_qc_flags(token, config)
    py_fn = _try_python_trim(path, token)
    dest_r1.parent.mkdir(parents=True, exist_ok=True)
    if py_fn is not None:
        result = py_fn(
            str(r1),
            str(r2) if r2 is not None and r2.exists() else None,
            str(dest_r1),
            str(dest_r2) if dest_r2 is not None else None,
            cfg,
        )
        if result:
            return [str(p) for p in result]
        return [str(dest_r1)] + ([str(dest_r2)] if dest_r2 is not None else [])
    cmd = [str(path), "-i", str(r1), "-o", str(dest_r1.parent)]
    if r2 is not None and r2.exists():
        cmd.extend(["-I", str(r2)])
    cmd.extend(list(cfg.get("extra_argv") or []))
    subprocess.check_call(cmd)
    return [str(dest_r1)] + ([str(dest_r2)] if dest_r2 is not None and dest_r2.exists() else [])


def trim_directory(
    src: Path,
    dest: Path,
    *,
    stage_qc: Optional[str] = None,
    postfix_map: Optional[Dict[str, str]] = None,
    config: Optional[Dict[str, Any]] = None,
) -> List[str]:
    """Trim every FASTQ pair in ``src`` into ``dest`` (identity when QC is unset)."""
    src_p = as_path(src)
    dest_p = as_path(dest)
    dest_p.mkdir(parents=True, exist_ok=True)
    cfg = dict(config or {})
    mapping = parse_qc_postfix(postfix_map if postfix_map is not None else cfg.get("qc_postfix"))
    written: List[str] = []
    for r1 in list_all_r1_files(src_p):
        sample = sample_name_from_r1(r1)
        r2 = r2_for_r1(r1)
        if not r2.exists():
            found = find_fastq_mate(src_p, sample, "R2")
            r2 = found
        dest_r1 = dest_p / r1.name
        dest_r2 = dest_p / r2.name if r2 is not None else None
        tool = qc_name_for_sample(sample, stage_qc=stage_qc, postfix_map=mapping)
        written.extend(
            run_qc_pair(tool, r1, r2, dest_r1, dest_r2, cfg)
        )
    return written


def trim_stage(output_dir: Path, stage: str, config: Optional[Dict[str, Any]] = None) -> List[str]:
    cfg = dict(config or {})
    root = as_path(output_dir)
    if stage == "generated":
        src = root / "regenerated"
        dest = root / "regenerated_trimmed"
        stage_qc = cfg.get("qc_generated")
        if stage_qc is None:
            stage_qc = cfg.get("qc")
    else:
        src = root / "initial"
        dest = root / "initial_trimmed"
        stage_qc = cfg.get("qc_initial")
        if stage_qc is None:
            stage_qc = cfg.get("qc")
    return trim_directory(
        src,
        dest,
        stage_qc=stage_qc,
        postfix_map=cfg.get("qc_postfix"),
        config=cfg,
    )


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="samovar.qc", description="Trim metagenome FASTQ")
    sub = parser.add_subparsers(dest="cmd", required=True)
    trim = sub.add_parser("trim", help="Trim a run stage or a directory")
    trim.add_argument("--src", default="", help="Source FASTQ directory")
    trim.add_argument("--dest", default="", help="Destination FASTQ directory")
    trim.add_argument("--output_dir", "--output-dir", dest="output_dir", default="")
    trim.add_argument("--stage", choices=("initial", "generated"), default="initial")
    trim.add_argument("--qc", default="", help="QC tool name (empty = identity)")
    trim.add_argument("--config", default="", help="YAML with qc / qc_postfix / flags")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    import yaml

    args = _parser().parse_args(list(argv) if argv is not None else None)
    cfg: Dict[str, Any] = {}
    if args.config:
        path = Path(args.config)
        if path.is_file():
            loaded = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
            if isinstance(loaded, dict):
                cfg.update(loaded)
    if args.qc:
        cfg["qc"] = args.qc
        if args.stage == "generated":
            cfg.setdefault("qc_generated", args.qc)
        else:
            cfg.setdefault("qc_initial", args.qc)
    if args.src and args.dest:
        trim_directory(
            Path(args.src),
            Path(args.dest),
            stage_qc=cfg.get("qc_generated" if args.stage == "generated" else "qc_initial", cfg.get("qc")),
            postfix_map=cfg.get("qc_postfix"),
            config=cfg,
        )
        return 0
    out = args.output_dir or cfg.get("output_dir")
    if not out:
        print("Error: --output_dir or --src/--dest is required", file=sys.stderr)
        return 2
    trim_stage(Path(out), args.stage, cfg)
    return 0


if __name__ == "__main__":
    sys.exit(main())
