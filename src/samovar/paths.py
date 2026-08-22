"""Install-layout helpers: repo root, XDG config, and executable lookup.

Config is read from (first hit wins):

1. ``$SAMOVAR_CONFIG``
2. ``$XDG_CONFIG_HOME/samovar/config.json`` (default ``~/.config/samovar/config.json``)
3. ``<repo>/build/config.json``
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

PACKAGE_VERSION = "0.10.10"

KNOWN_TOOLS = (
    "kraken2",
    "kaiju",
    "kraken",
    "krakenuniq",
    "metaphlan",
    "metaphlan4",
    "centrifuge",
    "iss",
    "opal.py",
    "opal",
    "R",
    "Rscript",
    "snakemake",
    "python",
    "python3",
)


def user_config_dir() -> Path:
    xdg = os.environ.get("XDG_CONFIG_HOME", "").strip()
    if xdg:
        return Path(xdg) / "samovar"
    return Path.home() / ".config" / "samovar"


def user_config_path() -> Path:
    override = os.environ.get("SAMOVAR_CONFIG", "").strip()
    if override:
        return Path(override)
    return user_config_dir() / "config.json"


def repo_root() -> Path:
    env = os.environ.get("SAMOVAR_ROOT", "").strip()
    if env:
        return Path(env).resolve()
    cfg = _read_json(user_config_path())
    root = (cfg or {}).get("root")
    if root and Path(root).is_dir():
        return Path(root).resolve()
    # src/samovar/paths.py → repo
    return Path(__file__).resolve().parent.parent.parent


def workflow_dir() -> Path:
    return repo_root() / "workflow"


def test_genomes_dir() -> Path:
    cfg = load_config()
    override = cfg.get("test_genomes") or os.environ.get("SAMOVAR_TEST_GENOMES", "")
    if override:
        return Path(override)
    return repo_root() / "data" / "test_genomes"


def absolute_path(path: Optional[str]) -> str:
    """Resolve ``path`` against the current working directory.

    Used at ``samovar prepare`` so generated YAML/scripts do not depend on cwd
    at ``samovar exec`` (e.g. Slurm ``cd $OUT``).
    """
    if path is None:
        return ""
    raw = str(path).strip()
    if not raw:
        return raw
    candidate = Path(raw).expanduser()
    if not candidate.is_absolute():
        candidate = Path.cwd() / candidate
    try:
        return str(candidate.resolve())
    except OSError:
        return str(candidate.absolute())


def _read_json(path: Path) -> Optional[Dict[str, Any]]:
    try:
        with open(path, encoding="utf-8") as handle:
            data = json.load(handle)
        return data if isinstance(data, dict) else None
    except (OSError, json.JSONDecodeError, TypeError):
        return None


def load_config() -> Dict[str, Any]:
    candidates = []
    override = os.environ.get("SAMOVAR_CONFIG", "").strip()
    if override:
        candidates.append(Path(override))
    candidates.append(user_config_dir() / "config.json")
    candidates.append(repo_root() / "build" / "config.json")
    merged: Dict[str, Any] = {}
    # Later files overlay earlier so repo build/ can fill gaps, user config wins.
    for path in reversed(list(dict.fromkeys(candidates))):
        data = _read_json(path)
        if data:
            merged.update(data)
    return merged


def python_path() -> str:
    cfg = load_config()
    configured = (cfg.get("python_path") or "").strip()
    if configured and Path(configured).is_file():
        return configured
    found = shutil.which(configured or "python3") or shutil.which("python")
    return found or sys.executable


def _split_path_value(value: Any) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple)):
        items: List[str] = []
        for part in value:
            items.extend(_split_path_value(part))
        return items
    text = str(value).strip()
    if not text or text.startswith("_"):
        return []
    return [
        piece.strip()
        for piece in text.replace(";", ":").split(":")
        if piece.strip() and piece.strip() not in {"$PATH", "${PATH}"}
    ]


def _dir_for_path_entry(raw: str, *, env_prefix: bool = False) -> Optional[str]:
    """Directory to prepend to PATH for an executable or a conda/module prefix."""
    text = str(raw).strip()
    if not text:
        return None
    path = Path(text).expanduser()
    # Bare command names live on PATH; they are not directories to prepend.
    if not path.is_absolute() and len(path.parts) == 1:
        return None
    if path.is_file():
        return str(path.resolve().parent)
    bindir = path / "bin"
    if path.is_dir():
        if path.name == "bin":
            return str(path.resolve())
        if bindir.is_dir():
            return str(bindir.resolve())
        return str(path.resolve())
    # Not on disk yet (other HPC / module not loaded). Keep .../bin as-is;
    # treat other env_prefix values as conda/module prefixes.
    if path.name == "bin":
        return str(path)
    if env_prefix:
        return str(bindir)
    if path.name in KNOWN_TOOLS or path.suffix:
        parent = path.parent
        return None if str(parent) in {".", ""} else str(parent)
    return str(path)


def collect_runtime_path_dirs(cfg: Optional[Dict[str, Any]] = None) -> List[str]:
    """Bin directories from config so bare ``bash .log/samovar.sh`` finds tools.

    Sources, in order:

    * repo ``bin/``
    * ``python_path`` / ``iss_path`` / ``r_path`` parent dirs
    * ``path`` / ``extra_path`` (string or list) — extra env ``bin/`` dirs
    * parent dirs of ``tools.*`` executables
    * ``tool_envs.<name>`` conda/module prefixes (``<prefix>/bin``)
    """
    cfg = dict(cfg or load_config())
    ordered: List[str] = []
    seen = set()

    def add(raw: Optional[str], env_prefix: bool = False) -> None:
        if not raw:
            return
        directory = _dir_for_path_entry(str(raw), env_prefix=env_prefix)
        if not directory or directory in seen:
            return
        seen.add(directory)
        ordered.append(directory)

    add(str(repo_root() / "bin"))
    add(cfg.get("python_path"))
    add(cfg.get("iss_path"))
    add(cfg.get("r_path"))
    for extra in _split_path_value(cfg.get("path") or cfg.get("extra_path")):
        add(extra, env_prefix=True)
    tools = cfg.get("tools") if isinstance(cfg.get("tools"), dict) else {}
    for key, value in tools.items():
        if str(key).startswith("_"):
            continue
        add(str(value) if value is not None else "")
    envs = cfg.get("tool_envs") if isinstance(cfg.get("tool_envs"), dict) else {}
    for key, value in envs.items():
        if str(key).startswith("_"):
            continue
        add(str(value) if value is not None else "", env_prefix=True)
    return ordered


def runtime_path_prefix(cfg: Optional[Dict[str, Any]] = None) -> str:
    return ":".join(collect_runtime_path_dirs(cfg))


def iss_executable() -> str:
    """ISS CLI: config ``iss_path``, then PATH."""
    cfg = load_config()
    configured = (cfg.get("iss_path") or "").strip() or "iss"
    resolved = resolve_executable(configured, tool_key="iss")
    token = (resolved or "iss").split()[0]
    if token and Path(token).is_file() and os.access(token, os.X_OK):
        return str(Path(token).resolve())
    found = shutil.which(token) or shutil.which("iss")
    return found or token or "iss"


def ncbi_email() -> str:
    for key in ("NCBI_EMAIL", "ENTREZ_EMAIL", "SAMOVAR_EMAIL"):
        value = os.environ.get(key, "").strip()
        if value:
            return value
    cfg = load_config()
    email = str(cfg.get("ncbi_email") or "").strip()
    if email:
        return email
    return "anonymous@example.com"


def resolve_executable(name_or_path: Optional[str], tool_key: Optional[str] = None) -> str:
    """Return an absolute executable path when possible; otherwise the original token."""
    if not name_or_path:
        return ""
    raw = str(name_or_path).strip()
    if not raw:
        return raw
    token = raw.split()[0]
    rest = raw[len(token) :].lstrip()

    cfg = load_config()
    tools = cfg.get("tools") if isinstance(cfg.get("tools"), dict) else {}
    name = Path(tool_key or token).name
    mapped = str(tools.get(name) or tools.get(Path(name).stem) or "").strip()
    envs = cfg.get("tool_envs") if isinstance(cfg.get("tool_envs"), dict) else {}
    env_root = str(envs.get(name) or envs.get(Path(name).stem) or "").strip()
    env_candidates: List[str] = []
    if env_root:
        root = Path(env_root).expanduser()
        env_candidates.append(str(root / "bin" / name))
        env_candidates.append(str(root / name))

    # Other-env prefixes in tool_envs beat a same-named binary already on PATH
    # (and beat tools.* discovered at install time).
    candidates = [*env_candidates, mapped, token]
    for cand in candidates:
        if not cand:
            continue
        path = Path(cand).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            resolved = str(path.resolve())
            return f"{resolved} {rest}".strip()
        if path.is_absolute():
            continue
        which = shutil.which(cand)
        if which:
            return f"{which} {rest}".strip()
    return raw


def write_config(data: Dict[str, Any], also_repo_build: bool = True) -> Path:
    """Write user config (and a copy under ``<repo>/build/config.json``)."""
    payload = dict(data)
    payload.setdefault("version", PACKAGE_VERSION)
    payload.setdefault("root", str(repo_root()))
    user_dir = user_config_dir()
    user_dir.mkdir(parents=True, exist_ok=True)
    dest = user_config_dir() / "config.json"
    dest.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    if also_repo_build:
        try:
            build = repo_root() / "build"
            build.mkdir(parents=True, exist_ok=True)
            (build / "config.json").write_text(
                dest.read_text(encoding="utf-8"), encoding="utf-8"
            )
        except OSError:
            # Read-only shared checkout: user config is enough.
            pass
    return dest


def discover_tools() -> Dict[str, str]:
    found: Dict[str, str] = {}
    for name in KNOWN_TOOLS:
        path = shutil.which(name)
        if path:
            found[name] = path
    return found


def annotation_regenerate_r() -> Optional[Path]:
    """Locate the optional R regenerator script.

    Lookup order (first hit wins):

    1. ``$SAMOVAR_R_REGENERATE`` / ``$SAMOVAR_ANNOTATION_REGENERATE_R``
    2. config key ``annotation_regenerate_r``
    3. ``~/.config/samovar/annotation_regenerate.R`` written by ``./install.sh R-package``

    There is no R driver under ``workflow/`` on the Python tree.
    """
    for key in ("SAMOVAR_R_REGENERATE", "SAMOVAR_ANNOTATION_REGENERATE_R"):
        env = os.environ.get(key, "").strip()
        if env:
            return Path(env)
    cfg = load_config()
    override = str(cfg.get("annotation_regenerate_r") or "").strip()
    if override:
        return Path(override)
    xdg = user_config_dir() / "annotation_regenerate.R"
    if xdg.is_file():
        return xdg
    return None


def cxx_compiler() -> Optional[str]:
    env = os.environ.get("CXX", "").strip()
    for name in (env, "g++", "c++", "clang++"):
        if not name:
            continue
        found = shutil.which(name)
        if found:
            return found
    return None


def smoke_test() -> List[str]:
    """Return a list of problem strings (empty means OK)."""
    problems: List[str] = []
    try:
        import samovar  # noqa: F401
    except Exception as exc:  # pragma: no cover
        problems.append(f"import samovar failed: {exc}")
    iss = iss_executable()
    if not (iss and (Path(iss).is_file() or shutil.which(iss) or shutil.which("iss"))):
        problems.append("iss CLI not on PATH (install InSilicoSeq / insilicoseq)")
    if not shutil.which("snakemake"):
        problems.append("snakemake not on PATH")
    root = repo_root()
    if not (root / "workflow" / "annotators" / "Snakefile").is_file():
        problems.append(f"workflow Snakefiles missing under {root}")
    if cxx_compiler() is None:
        problems.append("no C++ compiler (g++/c++/clang++) for annotation combiner")
    try:
        from samovar.viz_annotation import require_viz_backend

        require_viz_backend()
    except Exception as exc:
        problems.append(str(exc))
    return problems
