"""``samovar tools import`` — register a binary in the install config."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path
from typing import Optional

from samovar.main_config import (
    TOOL_GROUPS,
    iter_tools,
    normalize_tool_group,
    parse_tool_entry,
    set_tool,
)
from samovar.paths import load_config, update_config


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
    also_repo_build: bool = True,
) -> list:
    """Write ``tools.<name> = [env, workflow, path, group, flags?]`` and return the spec."""
    name = str(name or "").strip()
    if not name:
        raise ValueError("--name is required")
    group = normalize_tool_group(tool_type)
    env = str(env or "").strip()
    exe = str(exec_name or "").strip() or name
    path = resolve_import_path(name=name, exec_name=exe, exec_path=exec_path, env=env)
    workflow = env if env else "bash"
    cfg = load_config()
    set_tool(
        cfg,
        name,
        path=path,
        env=env,
        workflow=workflow,
        group=group,
        flags=str(flags or "").strip(),
    )
    spec = parse_tool_entry(iter_tools(cfg).get(name), name)
    update_config({"tools": {name: spec}}, also_repo_build=also_repo_build)
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
        help=f"Tool group: {', '.join(TOOL_GROUPS)} (aliases: a, reads, meta, table, score)",
    )
    parser.add_argument(
        "--flags",
        default="",
        help="Optional extra CLI flags stored as the 5th tools.* slot (omitted when empty)",
    )
    return parser


def main(argv: Optional[list] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        spec = import_tool(
            name=args.name,
            tool_type=args.type,
            env=args.env,
            exec_name=args.exec_name,
            exec_path=args.exec_path,
            flags=args.flags,
        )
    except (ValueError, FileNotFoundError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    env, workflow, path, group = spec[:4]
    extra = spec[4] if len(spec) > 4 else ""
    row = f"[{env!r}, {workflow!r}, {path}, {group}"
    if extra:
        row += f", {extra!r}"
    row += "]"
    print(f"Imported tools.{args.name} = {row}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
