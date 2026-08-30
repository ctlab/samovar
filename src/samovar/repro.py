"""Hydra run snapshots: capture → export (fold) → install/import (unfold) → reproduce.

Hydra (facebookresearch/hydra, hydra-core) composes YAML + CLI overrides and
saves ``.hydra/config.yaml`` + ``overrides.yaml``. It does **not** lock conda
packages or pack databases. This module adds the missing lockfile, relocatable
``lazy-install`` recipes, and ``samovar export`` / ``samovar reproduce``.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import subprocess
import sys
import tarfile
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from omegaconf import DictConfig, OmegaConf

from samovar.main_config import (
    TOOL_GROUP_BY_NAME,
    iter_tool_records,
    tool_group_for,
)
from samovar.paths import (
    _code_repo_root,
    add_output_dir_argument,
    argv_has_output_dir,
    load_config,
    sidecar_envs_dir,
    token_is_output_dir_flag,
    user_config_path,
)
from samovar.tool_spec import (
    bare_tool_name,
    join_tool_key,
    lazy_install_for,
    probe_tool_version,
    split_tool_key,
)
from samovar.version import get_version

HYDRA_DIRNAME = ".hydra"
CAPTURE_ENV = (
    "NCBI_EMAIL",
    "ENTREZ_EMAIL",
    "SAMOVAR_EMAIL",
    "SAMOVAR_ROOT",
    "SAMOVAR_CONFIG",
    "SAMOVAR_DATABASE",
    "SAMOVAR_TAXDUMP",
    "SAMOVAR_CONDA",
    "SAMOVAR_ENVS",
    "CONDA_PREFIX",
    "PYTHON_PATH",
    "PYTHONPATH",
    "PATH",
)
BUILTIN_NAMES = set(TOOL_GROUP_BY_NAME) | {
    "dummy",
    "constant9606",
    "dummy9606",
    "constant",
    "samovar",
}


def hydra_dir(output_dir: Path) -> Path:
    return Path(output_dir).expanduser().resolve() / HYDRA_DIRNAME


def snapshot_path(output_dir: Path) -> Path:
    return hydra_dir(output_dir) / "config.yaml"


def _oc_escape_script(text: str) -> str:
    """Keep bash ``${PREFIX}`` from being parsed as an OmegaConf resolver."""
    return str(text).replace("${", "$${")


def _plain(obj: Any, default: Any = None) -> Any:
    """OmegaConf.to_container that also accepts native list/dict (empty lists are falsy)."""
    if obj is None:
        return default
    if OmegaConf.is_config(obj):
        try:
            value = OmegaConf.to_container(obj, resolve=True)
        except Exception:
            value = OmegaConf.to_container(obj, resolve=False)
        return default if value is None else value
    return obj


def load_lazy_install_text(raw: str) -> str:
    """CLI ``--lazy-install``: command, ``@file.sh``, or ``-`` (stdin)."""
    text = str(raw or "")
    if not text.strip():
        return ""
    stripped = text.strip()
    if stripped == "-":
        return sys.stdin.read()
    path_hint = stripped[1:] if stripped.startswith("@") else stripped
    candidate = Path(path_hint).expanduser()
    look_like_file = stripped.startswith("@") or candidate.suffix in {
        ".sh",
        ".bash",
        ".txt",
    }
    if candidate.is_file() and (look_like_file or stripped.startswith("@")):
        return candidate.read_text(encoding="utf-8")
    if candidate.is_file() and "\n" not in stripped and " " not in stripped:
        return candidate.read_text(encoding="utf-8")
    return text


def is_builtin_tool(name: str) -> bool:
    return bare_tool_name(name) in BUILTIN_NAMES


def hydra_dir_from_user(path: Path) -> Path:
    """Accept a run dir, a ``.hydra`` dir, a config.yaml, or an extracted export."""
    raw = Path(path).expanduser()
    if raw.is_file() and raw.name in {"config.yaml", "samovar.yaml"}:
        return raw.parent
    if raw.is_dir() and (raw / "config.yaml").is_file():
        return raw
    if raw.is_dir() and (raw / HYDRA_DIRNAME / "config.yaml").is_file():
        return hydra_dir(raw)
    raise FileNotFoundError(f"No Hydra snapshot under {path}")


def _ns_to_dict(args: Any) -> Dict[str, Any]:
    if args is None:
        return {}
    if isinstance(args, dict):
        data = dict(args)
    else:
        data = dict(vars(args))
    out: Dict[str, Any] = {}
    for key, value in data.items():
        if str(key).startswith("_"):
            continue
        try:
            json.dumps(value)
        except TypeError:
            out[str(key)] = str(value)
            continue
        out[str(key)] = value
    return out


def _argv_list(argv: Optional[Sequence[str]] = None) -> List[str]:
    if argv is not None:
        return [str(x) for x in argv]
    return [str(x) for x in sys.argv[1:]]


def args_to_overrides(stage: str, args: Mapping[str, Any]) -> List[str]:
    """Hydra-style ``stage.key=value`` overrides (skip empty / noisy fields)."""
    skip = {"custom_cli_flags", "tool_flags"}
    rows: List[str] = []
    for key, value in args.items():
        if key in skip or value in (None, "", [], {}, False):
            continue
        if isinstance(value, (list, dict)):
            rows.append(f"{stage}.{key}={json.dumps(value)}")
        else:
            rows.append(f"{stage}.{key}={value}")
    return rows


def capture_env() -> Dict[str, Any]:
    from samovar.paths import python_path

    variables = {key: os.environ[key] for key in CAPTURE_ENV if os.environ.get(key)}
    return {
        "python": python_path(),
        "python_version": sys.version.replace("\n", " "),
        "hydra": _hydra_version(),
        "platform": platform.platform(),
        "machine": platform.machine(),
        "hostname": platform.node(),
        "cwd": str(Path.cwd()),
        "time": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "variables": variables,
        "conda_prefix": os.environ.get("CONDA_PREFIX", ""),
    }


def capture_tools(cfg: Optional[Mapping[str, Any]] = None) -> Tuple[Dict[str, Any], Dict[str, Any], List[str]]:
    data = cfg or load_config()
    tools: Dict[str, Any] = {}
    imports: Dict[str, Any] = {}
    warnings: List[str] = []
    for key, rec in iter_tool_records(data).items():
        if ":" not in str(key) and any(
            str(other).startswith(str(key) + ":") for other in iter_tool_records(data)
        ):
            continue
        bare, ver = split_tool_key(key)
        exe = dict(rec.get("exec") or {})
        path = str(exe.get("path") or "")
        found = ver or probe_tool_version(path, bare)
        if bare in {"samovar"}:
            found = get_version()
        lazy = str(rec.get("lazy-install") or "").strip() or lazy_install_for(bare, found)
        builtin = is_builtin_tool(bare)
        row = {
            "name": bare,
            "version": found,
            "type": rec.get("type") or tool_group_for(bare),
            "builtin": builtin,
            "exec": exe,
            "flags": rec.get("flags") or "",
            "flags_translate": dict(rec.get("flags-translate") or {}),
            "lazy_install": _oc_escape_script(lazy),
            "inputs": rec.get("inputs") or "",
        }
        tools[join_tool_key(bare, found) if found else bare] = row
        if not builtin:
            imports[bare] = row
            if not lazy:
                warnings.append(
                    f"custom tool {bare!r} has no lazy-install; export/reproduce cannot reinstall it"
                )
    return tools, imports, warnings


def capture_databases(cfg: Optional[Mapping[str, Any]] = None) -> Dict[str, List[Dict[str, str]]]:
    from samovar.db_spec import iter_database_records, official_url_for

    data = cfg or load_config()
    out: Dict[str, List[Dict[str, str]]] = {}
    for tool, grouped in iter_database_records(data).items():
        packed = []
        for rec in grouped.values():
            version = str(rec.get("_version") or rec.get("version") or "")
            packed.append(
                {
                    "name": str(rec.get("name") or ""),
                    "path": str(rec.get("path") or ""),
                    "flags": str(rec.get("flags") or ""),
                    "version": version,
                    "lazy_download": _oc_escape_script(str(rec.get("lazy-download") or "")),
                    "url": official_url_for(
                        tool, str(rec.get("name") or ""), version, str(rec.get("url") or "")
                    ),
                }
            )
        if packed:
            out[tool] = packed
    return out


def empty_snapshot() -> DictConfig:
    return OmegaConf.create(
        {
            "samovar_version": get_version(),
            "command": "",
            "stages": [],
            "generate": {},
            "prepare": {},
            "exec": {},
            "env": {},
            "tools": {},
            "imports": {},
            "databases": {},
            "genomes": {},
            "custom_flags": [],
            "warnings": [],
        }
    )


def compose_snapshot(hydra_root: Path, extra_overrides: Optional[Sequence[str]] = None) -> DictConfig:
    """Replay a run: overlay saved ``config.yaml``, then Hydra-style extra overrides.

    Hydra's job output is the composed snapshot; ``overrides.yaml`` is the audit
    trail of the original CLI. Reproduce applies *new* overrides (``key=value``)
    on top of that snapshot — it does not re-compose the original CLI.
    """
    saved = OmegaConf.load(hydra_root / "config.yaml")
    extra = [str(x) for x in (extra_overrides or []) if str(x).strip()]
    try:
        from hydra import compose, initialize_config_module
        from hydra.core.global_hydra import GlobalHydra

        if GlobalHydra.instance().is_initialized():
            GlobalHydra.instance().clear()
        with initialize_config_module(
            version_base="1.3",
            config_module="samovar.conf",
            job_name="samovar-reproduce",
        ):
            recomposed = compose(config_name="config", overrides=extra)
        merged = OmegaConf.merge(recomposed, saved)
    except Exception:
        merged = saved
    if extra:
        try:
            merged = OmegaConf.merge(merged, OmegaConf.from_dotlist(extra))
        except Exception:
            pass
    return merged


def record_stage(
    stage: str,
    output_dir: Optional[str],
    args: Any = None,
    argv: Optional[Sequence[str]] = None,
    extra: Optional[Mapping[str, Any]] = None,
) -> Optional[Path]:
    """Merge this generate/prepare/exec invocation into ``$out/.hydra``."""
    if not output_dir:
        return None
    out = Path(output_dir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    dest = hydra_dir(out)
    dest.mkdir(parents=True, exist_ok=True)
    snap_file = dest / "config.yaml"
    if snap_file.is_file():
        snap = OmegaConf.load(snap_file)
    else:
        snap = empty_snapshot()
    cfg = load_config()
    tools, imports, warnings = capture_tools(cfg)
    parsed = _ns_to_dict(args)
    argv_list = _argv_list(argv)
    if parsed and not argv_has_output_dir(argv_list):
        argv_list = _argv_from_args(stage, parsed)
    block = {"argv": argv_list, "args": parsed}
    OmegaConf.update(snap, stage, block, merge=True)
    stages = list(_plain(snap.get("stages"), []) or [])
    if stage not in stages:
        stages.append(stage)
    snap.stages = stages
    snap.command = stage
    snap.samovar_version = get_version()
    snap.env = capture_env()
    snap.tools = tools
    snap.imports = imports
    snap.databases = capture_databases(cfg)
    snap.custom_flags = list(cfg.get("custom-flags") or ["--threads", "--cores"])
    existing_warn = list(_plain(snap.get("warnings"), []) or [])
    for item in warnings:
        if item not in existing_warn:
            existing_warn.append(item)
    snap.warnings = existing_warn
    if extra:
        for key, value in extra.items():
            if value in (None, "", [], {}):
                continue
            OmegaConf.update(snap, str(key), value, merge=True)

    OmegaConf.save(snap, snap_file)
    overrides = args_to_overrides(stage, parsed)
    ov_file = dest / "overrides.yaml"
    previous = []
    if ov_file.is_file():
        previous = list(_plain(OmegaConf.load(ov_file), []) or [])
    OmegaConf.save(previous + overrides, ov_file)
    hydra_meta = {
        "hydra": {
            "run": {"dir": str(out.resolve())},
            "job": {"name": f"samovar-{stage}", "chdir": False},
            "runtime": {
                "cwd": str(Path.cwd()),
                "output_dir": str(out.resolve()),
                "directory": str(out.resolve()),
                "version": _hydra_version(),
            },
        }
    }
    OmegaConf.save(OmegaConf.create(hydra_meta), dest / "hydra.yaml")
    return dest


def _hydra_version() -> str:
    try:
        import hydra

        return str(getattr(hydra, "__version__", "") or "")
    except Exception:
        return ""


def compare_tools(snap: DictConfig) -> List[Dict[str, str]]:
    """Rows of {name, snapshot, installed, status, lazy_install}."""
    rows: List[Dict[str, str]] = []
    tools = _plain(snap.get("tools"), {}) or {}
    for key, rec in tools.items():
        if not isinstance(rec, dict):
            continue
        bare = str(rec.get("name") or bare_tool_name(key))
        want = str(rec.get("version") or split_tool_key(key)[1])
        path = str((rec.get("exec") or {}).get("path") or "")
        found = ""
        if bare in {"samovar"}:
            found = get_version()
        elif want:
            found = probe_tool_version(path, bare) if path else ""
            if not found:
                from samovar.main_config import lookup_tool_record

                live = lookup_tool_record(load_config(), bare)
                if live:
                    live_path = str((live.get("exec") or {}).get("path") or "")
                    found = probe_tool_version(live_path, bare)
                    found = found or str(live.get("_version") or "")
                    path = live_path or path
        status = "ok"
        if want and found and found != want and found.split(".")[:2] != want.split(".")[:2]:
            status = "mismatch"
        elif want and not found:
            status = "missing"
        ft = rec.get("flags_translate") or rec.get("flags-translate") or {}
        if isinstance(ft, dict):
            flags_translate = " ".join(f"{k}:{v}" for k, v in ft.items() if k and v)
        else:
            flags_translate = str(ft or "")
        rows.append(
            {
                "name": bare,
                "snapshot": want,
                "installed": found,
                "status": status,
                "lazy_install": str(rec.get("lazy_install") or ""),
                "builtin": "1" if rec.get("builtin") else "0",
                "path": path,
                "type": str(rec.get("type") or ""),
                "flags": str(rec.get("flags") or ""),
                "flags_translate": flags_translate,
            }
        )
    return rows


def _inject_conda_prefix(script: str, prefix: Path) -> str:
    text = script.strip()
    if text.startswith(("conda install", "mamba install", "micromamba install")):
        parts = text.split(None, 2)
        if len(parts) >= 2 and "-p" not in text and "--prefix" not in text:
            return f"{parts[0]} install -y -p {prefix} " + (parts[2] if len(parts) > 2 else "")
    return text


def ensure_tool_installed(row: Mapping[str, str], *, dry_run: bool = False) -> Optional[Path]:
    """Create a sidecar env and run lazy-install when versions mismatch."""
    if row.get("status") == "ok":
        return None
    lazy = str(row.get("lazy_install") or "").strip()
    name = row["name"]
    ver = row.get("snapshot") or "any"
    if not lazy:
        return None
    prefix = sidecar_envs_dir() / f"repro-{name}-{ver}"
    script = _inject_conda_prefix(lazy, prefix)
    if dry_run:
        print(f"[dry-run] {name}: mkdir {prefix}")
        print(f"[dry-run] bash <<'EOF'\n{script}\nEOF")
        return prefix
    prefix.mkdir(parents=True, exist_ok=True)
    script_path = prefix / "lazy-install.sh"
    script_path.write_text(script + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["SAMOVAR_REPRO_PREFIX"] = str(prefix)
    env["PREFIX"] = str(prefix)
    print(f"Installing {name} {ver} into {prefix}")
    subprocess.check_call(["bash", str(script_path)], cwd=str(prefix), env=env)
    return prefix


def reimport_from_prefix(row: Mapping[str, str], prefix: Path) -> None:
    from samovar.tools_import import import_tool

    name = row["name"]
    binary = prefix / "bin" / name
    path = str(binary) if binary.is_file() else str(prefix)
    import_tool(
        name=name,
        tool_type=str(row.get("type") or "runtime"),
        env="conda" if (prefix / "conda-meta").is_dir() else "",
        exec_path=path,
        flags=str(row.get("flags") or ""),
        flags_translate=str(row.get("flags_translate") or ""),
        lazy_install=str(row.get("lazy_install") or ""),
        version=str(row.get("snapshot") or ""),
    )


def patch_argv_overrides(argv: Sequence[str], stage: str, extra: Sequence[str]) -> List[str]:
    """Apply Hydra ``stage.key=value`` overrides onto a recorded CLI argv."""
    updates: Dict[str, str] = {}
    prefix = f"{stage}."
    for item in extra:
        if "=" not in str(item):
            continue
        key, val = str(item).split("=", 1)
        if not key.startswith(prefix):
            continue
        key = key[len(prefix) :]
        if key.startswith("args."):
            key = key[5:]
        if key in {"argv", "args"}:
            continue
        updates[key] = val
    out = list(argv)
    for key, val in updates.items():
        flag = f"--{key}" if not str(key).startswith("-") else str(key)
        if flag in out:
            idx = out.index(flag)
            if idx + 1 < len(out) and not str(out[idx + 1]).startswith("-"):
                out[idx + 1] = str(val)
            else:
                out.insert(idx + 1, str(val))
        else:
            out.extend([flag, str(val)])
    return out


def rewrite_output_dir_argv(argv: Sequence[str], output_dir: str) -> List[str]:
    out: List[str] = []
    i = 0
    replaced = False
    tokens = list(argv)
    while i < len(tokens):
        tok = tokens[i]
        if "=" in tok and token_is_output_dir_flag(tok):
            out.extend(["--output_dir", output_dir])
            i += 1
            replaced = True
            continue
        if token_is_output_dir_flag(tok):
            out.extend(["--output_dir", output_dir])
            i += 2
            replaced = True
            continue
        out.append(tok)
        i += 1
    if not replaced:
        out.extend(["--output_dir", output_dir])
    return out


def run_recorded_stage(
    stage: str,
    block: Mapping[str, Any],
    output_dir: str,
    *,
    dry_run: bool = False,
) -> int:
    argv = list(block.get("argv") or [])
    if not argv:
        args = dict(block.get("args") or {})
        argv = _argv_from_args(stage, args)
    argv = rewrite_output_dir_argv(argv, output_dir)
    used = list(block.get("used_genomes") or [])
    if stage == "generate" and used:
        tokens = []
        for rec in used:
            if isinstance(rec, dict):
                if rec.get("host"):
                    continue
                tokens.append(str(rec.get("taxid") or rec.get("file") or rec.get("path") or ""))
            else:
                tokens.append(str(rec))
        tokens = [t for t in tokens if t]
        if tokens and "--genomes" not in argv:
            for tok in tokens:
                argv.extend(["--genomes", tok])
    root = _reproduce_code_root()
    samovar = root / "bin" / "samovar"
    if not samovar.is_file():
        which = shutil.which("samovar")
        samovar = Path(which) if which else samovar
    cmd = [str(samovar), stage, *argv]
    if stage == "prepare":
        cmd = [str(samovar), "prepare", *argv]
    print(" ".join(cmd))
    if dry_run:
        return 0
    env = os.environ.copy()
    env["PYTHONPATH"] = str(root / "src") + (
        os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else ""
    )
    env["PATH"] = str(root / "bin") + os.pathsep + env.get("PATH", "")
    proc = subprocess.run(cmd, env=env)
    return int(proc.returncode)


def _reproduce_code_root() -> Path:
    """Checkout that actually ships ``bin/samovar`` (not a config.json data root)."""
    env = os.environ.get("SAMOVAR_ROOT", "").strip()
    if env:
        root = Path(env).expanduser()
        if (root / "bin" / "samovar").is_file():
            return root.resolve()
    return _code_repo_root()


def _argv_from_args(stage: str, args: Mapping[str, Any]) -> List[str]:
    out: List[str] = []
    for key, value in args.items():
        if value in (None, "", False, [], {}) or key in {"custom_cli_flags", "tool_flags"}:
            continue
        flag = f"--{key}" if not str(key).startswith("-") else str(key)
        if value is True:
            out.append(flag)
            continue
        if isinstance(value, list):
            if value and isinstance(value[0], list):
                out.extend([flag, " ".join(str(x) for x in value[0])])
            else:
                for item in value:
                    out.extend([flag, str(item)])
            continue
        out.extend([flag, str(value)])
    return out


def database_src_exists(src: Path) -> bool:
    """True for a directory, file, or centrifuge-style prefix (``stem.1.cf``)."""
    if src.exists():
        return True
    return Path(str(src) + ".1.cf").is_file()


def copy_database_tree(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dest, dirs_exist_ok=True)
        return
    if src.is_file():
        shutil.copy2(src, dest)
        return
    parent = src.parent
    stem = src.name
    if not parent.is_dir():
        return
    dest.mkdir(parents=True, exist_ok=True)
    for child in parent.iterdir():
        if child.name == stem or child.name.startswith(stem + "."):
            shutil.copy2(child, dest / child.name)


def export_run(
    source: Path,
    archive: Path,
    *,
    mode: str = "full",
    include_genomes: bool = False,
) -> Path:
    hydra_root = hydra_dir_from_user(source)
    snap = OmegaConf.load(hydra_root / "config.yaml")
    warnings = list(_plain(snap.get("warnings"), []) or [])
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp) / "samovar-run"
        (root / HYDRA_DIRNAME).mkdir(parents=True)
        for name in ("config.yaml", "overrides.yaml", "hydra.yaml"):
            src = hydra_root / name
            if src.is_file():
                shutil.copy2(src, root / HYDRA_DIRNAME / name)
        payload = root / "payload"
        payload.mkdir()
        db_map: List[str] = []
        if mode == "full":
            databases = _plain(snap.get("databases"), {}) or {}
            for tool, rows in databases.items():
                for row in rows or []:
                    path = Path(str(row.get("path") or "")).expanduser()
                    name = str(row.get("name") or path.name or "db")
                    if not database_src_exists(path) and not str(row.get("lazy_download") or row.get("lazy-download") or "").strip():
                        warnings.append(f"database {tool}/{name} missing on disk: {path}")
                        continue
                    rel = Path("databases") / tool / name
                    db_map.append(f"{tool}\t{name}\t{rel.as_posix()}\t{row.get('flags') or ''}")
                    lazy = str(row.get("lazy_download") or row.get("lazy-download") or "").strip()
                    meta_dir = payload / "db-meta" / tool / name
                    meta_dir.mkdir(parents=True, exist_ok=True)
                    meta = {
                        "tool": tool,
                        "name": name,
                        "version": str(row.get("version") or ""),
                        "flags": row.get("flags") or "",
                        "url": str(row.get("url") or ""),
                    }
                    (meta_dir / "meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
                    if lazy:
                        (meta_dir / "lazy-download.sh").write_text(lazy.rstrip() + "\n", encoding="utf-8")
                    else:
                        warnings.append(
                            f"database {tool}/{name} has no lazy-download; "
                            "export does not pack database files — reproduce cannot re-fetch it"
                        )
            imports = _plain(snap.get("imports"), {}) or {}
            for name, rec in imports.items():
                rec = rec or {}
                lazy = str(rec.get("lazy_install") or "").strip()
                dest = payload / "tools" / name
                dest.mkdir(parents=True, exist_ok=True)
                meta = {
                    "name": name,
                    "type": rec.get("type") or "annotator",
                    "version": rec.get("version") or "",
                    "flags": rec.get("flags") or "",
                    "flags_translate": rec.get("flags_translate") or rec.get("flags-translate") or {},
                    "builtin": bool(rec.get("builtin")),
                }
                (dest / "meta.json").write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
                if not lazy:
                    warnings.append(
                        f"custom tool {name!r} has no lazy-install; "
                        "reproduce will not be able to rebuild it"
                    )
                    continue
                (dest / "lazy-install.sh").write_text(lazy.rstrip() + "\n", encoding="utf-8")
            if include_genomes:
                warnings.append(
                    "genome FASTA libraries are not packed; Hydra lists used taxids "
                    "and reproduce re-resolves them via --genomes / lazy catalog fetch"
                )
        (root / "warnings.txt").write_text("\n".join(warnings) + ("\n" if warnings else ""), encoding="utf-8")
        (payload / "databases.tsv").write_text("\n".join(db_map) + ("\n" if db_map else ""), encoding="utf-8")
        (root / "install.sh").write_text(_bundle_install_script(), encoding="utf-8")
        os.chmod(root / "install.sh", 0o755)
        (root / "README.txt").write_text(_bundle_readme(mode), encoding="utf-8")
        archive = Path(archive).expanduser()
        archive.parent.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "w:gz") as tar:
            tar.add(root, arcname="samovar-run")
    return archive


def _bundle_install_script() -> str:
    return """#!/bin/bash
# Unfold a SamovaR export: lazy-install custom tools, re-register databases, import.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SAMOVAR="${SAMOVAR:-samovar}"
export PYTHONPATH="${SAMOVAR_ROOT:+$SAMOVAR_ROOT/src:}${PYTHONPATH:-}"
if [ ! -x "$(command -v "$SAMOVAR" || true)" ]; then
    echo "samovar CLI not found. Install SamovaR first (./install.sh), then re-run."
    exit 1
fi
if [ -d "$ROOT/payload/tools" ]; then
    for dir in "$ROOT/payload/tools"/*; do
        [ -d "$dir" ] || continue
        name="$(basename "$dir")"
        echo "lazy-install $name"
        PREFIX="${SAMOVAR_ENVS:-$ROOT/envs}/$name"
        mkdir -p "$PREFIX"
        export PREFIX SAMOVAR_REPRO_PREFIX="$PREFIX"
        if [ -f "$dir/lazy-install.sh" ]; then
            bash "$dir/lazy-install.sh"
        else
            echo "WARNING: no lazy-install.sh for $name (custom tool cannot be rebuilt)"
        fi
        typ="annotator"
        ver=""
        if [ -f "$dir/meta.json" ]; then
            typ="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("type") or "annotator")' "$dir/meta.json")"
            ver="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("version") or "")' "$dir/meta.json")"
        fi
        if [ -x "$PREFIX/bin/$name" ]; then
            extra=()
            [ -n "$ver" ] && extra+=(--tool-version "$ver")
            "$SAMOVAR" tools import -n "$name" --type "$typ" --exec-path "$PREFIX/bin/$name" \\
                --lazy-install-file "$dir/lazy-install.sh" "${extra[@]}" || true
        fi
    done
fi
if [ -f "$ROOT/payload/databases.tsv" ]; then
    while IFS=$'\\t' read -r tool name rel flags; do
        [ -n "$tool" ] || continue
        dest="$ROOT/payload/$rel"
        extra=()
        [ -n "${flags:-}" ] && extra+=(--flags "$flags")
        if [ -f "$ROOT/payload/db-meta/$tool/$name/lazy-download.sh" ]; then
            extra+=(--lazy-download-file "$ROOT/payload/db-meta/$tool/$name/lazy-download.sh")
            if [ ! -e "$dest" ]; then
                mkdir -p "$dest"
                PREFIX="$dest" bash "$ROOT/payload/db-meta/$tool/$name/lazy-download.sh" || true
            fi
        fi
        if [ -f "$ROOT/payload/db-meta/$tool/$name/meta.json" ]; then
            ver="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("version") or "")' "$ROOT/payload/db-meta/$tool/$name/meta.json")"
            url="$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1])).get("url") or "")' "$ROOT/payload/db-meta/$tool/$name/meta.json")"
            [ -n "$ver" ] && extra+=(--tool-version "$ver")
            [ -n "$url" ] && extra+=(--url "$url")
        fi
        "$SAMOVAR" tools import -n "$name" --type database --tool "$tool" --exec-path "$dest" "${extra[@]}" || true
    done < "$ROOT/payload/databases.tsv"
fi
echo "Unfold complete. Reproduce with:"
echo "  samovar reproduce $ROOT --output_dir /path/to/new_run"
"""


def _bundle_readme(mode: str) -> str:
    return (
        "SamovaR Hydra run bundle\n"
        f"mode: {mode}\n\n"
        "Hydra saved the composed config in .hydra/ (config.yaml, overrides.yaml).\n"
        "Used genomes are listed under genomes.used (not packed as FASTA).\n"
        "Databases are recorded as lazy-download recipes, not copied.\n\n"
        "  bash install.sh\n"
        "  samovar reproduce . --output_dir ./rerun\n"
    )


def extract_archive(archive: Path, dest: Path) -> Path:
    dest.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive, "r:*") as tar:
        kwargs = {}
        if hasattr(tarfile, "data_filter"):
            kwargs["filter"] = "data"
        tar.extractall(dest, **kwargs)
    inner = dest / "samovar-run"
    if (inner / HYDRA_DIRNAME / "config.yaml").is_file():
        return inner
    if (dest / HYDRA_DIRNAME / "config.yaml").is_file():
        return dest
    raise FileNotFoundError(f"archive {archive} has no .hydra/config.yaml")


def reproduce(
    source: Path,
    output_dir: str,
    *,
    extra_overrides: Optional[Sequence[str]] = None,
    install: bool = True,
    dry_run: bool = False,
    stages: Optional[Sequence[str]] = None,
) -> int:
    path = Path(source).expanduser()
    work = path
    tmp: Optional[tempfile.TemporaryDirectory] = None
    if path.is_file() and str(path).endswith((".tgz", ".tar.gz", ".tar")):
        tmp = tempfile.TemporaryDirectory()
        work = extract_archive(path, Path(tmp.name))
    hydra_root = hydra_dir_from_user(work)
    snap = compose_snapshot(hydra_root, extra_overrides)
    rows = compare_tools(snap)
    print("Tool versions vs this machine:")
    for row in rows:
        mark = row["status"]
        print(f"  [{mark:8}] {row['name']}: snapshot={row['snapshot'] or '?'} installed={row['installed'] or '—'}")
        if mark != "ok" and not row.get("lazy_install") and row.get("builtin") != "1":
            print(f"    WARNING: no lazy-install for custom tool {row['name']}")
    if install:
        for row in rows:
            if row["status"] == "ok":
                continue
            if not row.get("lazy_install"):
                if row.get("builtin") == "1":
                    print(f"WARNING: builtin {row['name']} version differs; continuing (no sidecar recipe).")
                    continue
                print(f"ERROR: cannot reinstall custom {row['name']} without lazy-install", file=sys.stderr)
                if tmp:
                    tmp.cleanup()
                return 1
            prefix = ensure_tool_installed(row, dry_run=dry_run)
            if prefix is not None and not dry_run:
                try:
                    reimport_from_prefix(row, prefix)
                except Exception as exc:
                    print(f"Warning: import after install failed for {row['name']}: {exc}", file=sys.stderr)
    wanted = list(stages or _plain(snap.get("stages"), []) or [])
    if not wanted:
        wanted = [str(snap.get("command") or "")]
    wanted = [s for s in wanted if s]
    extra = [str(x) for x in (extra_overrides or []) if str(x).strip()]
    code = 0
    for stage in wanted:
        block = _plain(snap.get(stage), {}) or {}
        if not block:
            continue
        args = dict(block.get("args") or {})
        for key, value in block.items():
            if key not in {"argv", "args"}:
                args[key] = value
        block = dict(block)
        block["args"] = args
        if extra:
            argv0 = list(block.get("argv") or _argv_from_args(stage, args))
            block["argv"] = patch_argv_overrides(argv0, stage, extra)
        genomes = _plain(snap.get("genomes"), {}) or {}
        if stage == "generate" and isinstance(genomes, dict):
            block["used_genomes"] = genomes.get("used") or []
        code = run_recorded_stage(stage, block, output_dir, dry_run=dry_run)
        if code != 0:
            break
    if tmp:
        tmp.cleanup()
    return code


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="samovar.repro")
    sub = parser.add_subparsers(dest="cmd", required=True)
    rec = sub.add_parser("record", help="Write .hydra snapshot for a stage")
    rec.add_argument("stage", choices=("generate", "prepare", "exec"))
    add_output_dir_argument(rec, required=True)
    rec.add_argument("argv", nargs=argparse.REMAINDER, help="Original CLI tokens after --")
    exp = sub.add_parser("export", help="Fold a run into a relocatable archive")
    exp.add_argument("source", nargs="?", default=".", help="Run dir, .hydra dir, or config.yaml")
    exp.add_argument("-o", "--output", default="samovar-run.tgz")
    exp.add_argument("--config", choices=("full", "slim"), default="full")
    exp.add_argument("--include-genomes", action="store_true")
    rep = sub.add_parser("reproduce", help="Unfold snapshot/archive and re-run")
    rep.add_argument("source", help="Run dir, .hydra, config.yaml, or .tgz")
    add_output_dir_argument(rep, required=True)
    rep.add_argument("--dry-run", action="store_true")
    rep.add_argument("--no-install", action="store_true")
    rep.add_argument(
        "--stages",
        default=None,
        help="Comma-separated stages to replay (generate,prepare,exec). Default: all recorded.",
    )
    rep.add_argument(
        "overrides",
        nargs="*",
        help="Hydra overrides, e.g. generate.n_samples=2",
    )
    return parser


def _stages_and_overrides(
    stages_arg: Optional[str],
    overrides: Optional[Sequence[str]],
) -> Tuple[Optional[List[str]], List[str]]:
    extra = [str(x) for x in (overrides or []) if str(x).strip()]
    wanted: List[str] = []
    if stages_arg:
        for part in str(stages_arg).split(","):
            part = part.strip()
            if not part:
                continue
            if "=" in part:
                extra.append(part)
            else:
                wanted.append(part)
    return (wanted or None), extra


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(list(argv) if argv is not None else None)
    if args.cmd == "record":
        argv = list(args.argv or [])
        if argv and argv[0] == "--":
            argv = argv[1:]
        record_stage(args.stage, args.output_dir, args=None, argv=argv)
        return 0
    if args.cmd == "export":
        path = export_run(
            Path(args.source),
            Path(args.output),
            mode=args.config,
            include_genomes=args.include_genomes,
        )
        print(f"Wrote {path}")
        return 0
    if args.cmd == "reproduce":
        stages, extra = _stages_and_overrides(args.stages, args.overrides)
        return reproduce(
            Path(args.source),
            args.output_dir,
            extra_overrides=extra,
            install=not args.no_install,
            dry_run=args.dry_run,
            stages=stages,
        )
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
