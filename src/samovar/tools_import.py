"""``samovar tools import`` — register a binary in the install config."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path
from typing import Optional

from samovar.main_config import (
    DEFAULT_SCORING_INPUTS,
    TOOL_GROUPS,
    iter_tools,
    lookup_tool_record,
    normalize_tool_group,
    parse_tool_entry,
    set_tool,
)
from samovar.paths import load_config, update_config
from samovar.tool_contracts import format_contract, run_contract_pytest


def resolve_import_path(
    *,
    name: str,
    exec_name: str,
    exec_path: str,
    env: str,
) -> str:
    """Return the path stored in ``tools.<name>``."""
    raw = str(exec_path or "").strip()
    exe = str(exec_name or name).strip() or name
    if raw:
        candidate = Path(raw).expanduser()
        if candidate.is_file():
            return str(candidate.resolve())
        if candidate.is_dir():
            for nested in (candidate / "bin" / exe, candidate / exe):
                try:
                    if nested.is_file():
                        return str(nested.resolve()) if env.lower() != "conda" else str(candidate.resolve())
                except OSError:
                    continue
            if env.lower() == "conda":
                return str(candidate.resolve())
            raise FileNotFoundError(
                f"--exec-path {raw} is a directory but {exe} was not found under bin/"
            )
        if candidate.is_absolute():
            raise FileNotFoundError(f"--exec-path not found: {raw}")
        which = shutil.which(raw)
        if which:
            return str(Path(which).resolve())
        raise FileNotFoundError(f"--exec-path not found on PATH: {raw}")
    which = shutil.which(exe) or shutil.which(name)
    if which:
        return str(Path(which).resolve())
    raise FileNotFoundError(
        f"Could not resolve executable {exe!r}. Pass --exec-path or put it on PATH."
    )


def import_tool(
    *,
    name: str,
    tool_type: str,
    env: str = "",
    exec_name: str = "",
    exec_path: str = "",
    flags: str = "",
    inputs: str = "",
    lazy_install: str = "",
    flags_translate: str = "",
    version: str = "",
    also_repo_build: bool = True,
) -> list:
    """Write ``tools.<name:version>`` object record; return list spec for callers."""
    name = str(name or "").strip()
    if not name:
        raise ValueError("--name is required")
    group = normalize_tool_group(tool_type)
    env = str(env or "").strip()
    exe = str(exec_name or "").strip() or name
    path = resolve_import_path(name=name, exec_name=exe, exec_path=exec_path, env=env)
    workflow = env if env else "bash"
    glob = str(inputs or "").strip()
    if group == "scoring" and not glob:
        glob = DEFAULT_SCORING_INPUTS
    from samovar.repro import load_lazy_install_text

    lazy = load_lazy_install_text(str(lazy_install or ""))
    cfg = load_config()
    kwargs = dict(
        path=path,
        env=env,
        workflow=workflow,
        group=group,
        flags=str(flags or "").strip(),
        lazy_install=lazy or None,
        flags_translate=flags_translate or None,
        version=str(version or "").strip() or None,
    )
    if glob:
        kwargs["inputs"] = glob
    set_tool(cfg, name, **kwargs)
    from samovar.tool_spec import bare_tool_name

    stored = dict.get(cfg, "tools") or {}
    disk_key = name
    for key in stored:
        if bare_tool_name(key) == name:
            disk_key = key
            break
    rec = stored.get(disk_key)
    if not isinstance(rec, dict):
        rec = lookup_tool_record(cfg, name) or {}
    update_config({"tools": {disk_key: rec}}, also_repo_build=also_repo_build)
    return parse_tool_entry(iter_tools(load_config()).get(name), name)


DATABASE_TYPES = {"database", "db", "databases"}


def is_database_type(value: str) -> bool:
    return str(value or "").strip().lower().replace("-", "_") in DATABASE_TYPES


def resolve_database_path(exec_path: str) -> str:
    raw = str(exec_path or "").strip()
    if not raw:
        raise FileNotFoundError("--exec-path is required for --type database")
    candidate = Path(raw).expanduser()
    if candidate.exists():
        return str(candidate.resolve())
    prefix = Path(str(candidate) + ".1.cf")
    if prefix.is_file():
        return str(candidate)
    parent = candidate.parent
    stem = candidate.name
    if parent.is_dir():
        try:
            for child in parent.iterdir():
                if child.name == stem or child.name.startswith(stem + "."):
                    return str(candidate)
        except OSError:
            pass
    raise FileNotFoundError(f"--exec-path not found: {raw}")


def import_database(
    *,
    name: str,
    tool: str = "",
    exec_path: str = "",
    flags: str = "",
    lazy_download: str = "",
    version: str = "",
    url: str = "",
    also_repo_build: bool = True,
) -> dict:
    """Write ``databases.<tool>.<name:version>``; return the object record."""
    from samovar.db_spec import lookup_database_record, parse_tool_and_name
    from samovar.main_config import set_database
    from samovar.repro import load_lazy_install_text

    annotator, db_name = parse_tool_and_name(name, tool)
    if not annotator:
        raise ValueError("--tool is required for --type database (or pass -n kraken2:dbname)")
    if not db_name:
        raise ValueError("--name is required")
    path = resolve_database_path(exec_path)
    lazy = load_lazy_install_text(str(lazy_download or ""))
    cfg = load_config()
    rec = set_database(
        cfg,
        annotator,
        db_name,
        path=path,
        flags=str(flags or "").strip(),
        lazy_download=lazy or None,
        version=str(version or "").strip() or None,
        url=str(url or "").strip() or None,
    )
    update_config({"databases": cfg.get("databases") or {}}, also_repo_build=also_repo_build)
    stored = lookup_database_record(load_config(), annotator, rec.get("name") or db_name) or rec
    return stored


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="samovar tools import",
        description="Register an executable in the SamovaR install config (tools.*).",
    )
    parser.add_argument(
        "-n",
        "--name",
        required=True,
        help="Config key (e.g. kaiju). Used as --<name>-test in prepare.",
    )
    parser.add_argument(
        "--env",
        default="",
        help='Runtime env: empty (host PATH) or "conda" (sidecar prefix)',
    )
    parser.add_argument(
        "--exec",
        dest="exec_name",
        default="",
        help="Binary name (default: --name). Used inside a conda prefix bin/.",
    )
    parser.add_argument(
        "--exec-path",
        dest="exec_path",
        default="",
        help="File, conda prefix, or name on PATH. Default: resolve --exec / --name.",
    )
    parser.add_argument(
        "-t",
        "--type",
        required=True,
        help=(
            f"Tool group: {', '.join(TOOL_GROUPS)} "
            "(aliases: a, reads, meta, table, table-scoring, score, viz, ml); "
            "or 'database' / 'db' to register an annotator index under databases.*"
        ),
    )
    parser.add_argument(
        "--tool",
        default="",
        help="Annotator that uses this DB (required with --type database). "
        "Also accepted as -n tool:dbname.",
    )
    parser.add_argument(
        "--url",
        default="",
        help="Official archive URL stored on the database record (feeds lazy-download).",
    )
    parser.add_argument(
        "--flags",
        default="",
        help=(
            "Native CLI flags stored on the record. For --type database they are "
            "merged into the annotator's extra at prepare (e.g. --memory-mapping)."
        )
    )
    parser.add_argument(
        "--flags-translate",
        dest="flags_translate",
        default="",
        help='Map SamovaR flags to the tool CLI, e.g. "--threads:--threads --cores:-t"',
    )
    parser.add_argument(
        "--lazy-install",
        dest="lazy_install",
        default="",
        help=(
            'Install recipe for reproduce/export: a command, a @path to a bash file, '
            'or "-" for stdin. Multiline scripts are stored on the tool record.'
        ),
    )
    parser.add_argument(
        "--lazy-install-file",
        dest="lazy_install_file",
        default="",
        help="Read a bash file into lazy-install (same as --lazy-install @FILE).",
    )
    parser.add_argument(
        "--lazy-download",
        dest="lazy_download",
        default="",
        help="Database rebuild recipe (same syntax as --lazy-install). Used with --type database.",
    )
    parser.add_argument(
        "--lazy-download-file",
        dest="lazy_download_file",
        default="",
        help="Read a bash file into lazy-download (same as --lazy-download @FILE).",
    )
    parser.add_argument(
        "--tool-version",
        dest="version",
        default="",
        help="Version stored in the tools key (name:version). Default: probe --version.",
    )
    parser.add_argument(
        "--inputs",
        "--input",
        "--glob",
        dest="inputs",
        default="",
        help=(
            "Input glob under the run output dir (6th tools.* slot). "
            "Scoring/viz (--type scoring) default: *annotations. "
            "table-scoring does not use this glob. "
            "Examples: *annotations, *annotations/combined_annotation_table.csv, *_plots"
        ),
    )
    parser.add_argument(
        "--pytest",
        action="store_true",
        help=(
            "Run the in→out contract pytest for --type on --exec-path before "
            "writing the config. Import only if the test passes."
        ),
    )
    return parser


def main(argv: Optional[list] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if is_database_type(args.type):
            lazy = args.lazy_download or args.lazy_install
            extra_file = str(
                getattr(args, "lazy_download_file", "")
                or getattr(args, "lazy_install_file", "")
                or ""
            ).strip()
            if extra_file:
                lazy = "@" + extra_file
            rec = import_database(
                name=args.name,
                tool=str(getattr(args, "tool", "") or ""),
                exec_path=args.exec_path,
                flags=args.flags,
                lazy_download=lazy,
                version=args.version,
                url=str(getattr(args, "url", "") or ""),
            )
            print(f"Imported databases.{rec.get('tool')}.{rec.get('name')} = {rec}")
            return 0
        dest = resolve_import_path(
            name=args.name,
            exec_name=str(args.exec_name or "").strip() or args.name,
            exec_path=args.exec_path,
            env=str(args.env or ""),
        )
        group = normalize_tool_group(args.type)
        if args.pytest:
            print(f"Running contract pytest for {group} on {dest} …")
            try:
                code, output = run_contract_pytest(dest, group)
            except ValueError as exc:
                print(f"Error: {exc}", file=sys.stderr)
                print(format_contract(group), file=sys.stderr)
                return 1
            if code != 0:
                print(output, file=sys.stderr)
                print("Import blocked: contract pytest failed.", file=sys.stderr)
                print(format_contract(group), file=sys.stderr)
                return 1
            print(output.rstrip() or "contract pytest passed")
        lazy = args.lazy_install
        extra_file = str(getattr(args, "lazy_install_file", "") or "").strip()
        if extra_file:
            lazy = "@" + extra_file
        spec = import_tool(
            name=args.name,
            tool_type=args.type,
            env=args.env,
            exec_name=args.exec_name,
            exec_path=dest,
            flags=args.flags,
            inputs=args.inputs,
            lazy_install=lazy,
            flags_translate=args.flags_translate,
            version=args.version,
        )
    except (ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    from samovar.main_config import lookup_tool_record
    from samovar.paths import load_config

    rec = lookup_tool_record(load_config(), args.name) or {}
    print(f"Imported tools.{args.name} = {rec}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
