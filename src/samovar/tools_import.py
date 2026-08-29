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
    cfg = load_config()
    kwargs = dict(
        path=path,
        env=env,
        workflow=workflow,
        group=group,
        flags=str(flags or "").strip(),
        lazy_install=str(lazy_install or "").strip() or None,
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
            "(aliases: a, reads, meta, table, table-scoring, score, viz, ml)"
        ),
    )
    parser.add_argument(
        "--flags",
        default="",
        help="Optional extra CLI flags stored on the tool record",
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
        help='Install recipe stored on the tool, e.g. "conda install bioconda::kraken2"',
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
        spec = import_tool(
            name=args.name,
            tool_type=args.type,
            env=args.env,
            exec_name=args.exec_name,
            exec_path=dest,
            flags=args.flags,
            inputs=args.inputs,
            lazy_install=args.lazy_install,
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
