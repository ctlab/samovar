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

PACKAGE_VERSION = "0.10.4"

KNOWN_TOOLS = (
    "kraken2",
    "kaiju",
    "kraken",
    "krakenuniq",
    "metaphlan",
    "metaphlan4",
    "centrifuge",
    "iss",
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
    key = tool_key or Path(token).name
    mapped = str(tools.get(key) or tools.get(Path(token).stem) or "").strip()

    candidates = [token, mapped]
    for cand in candidates:
        if not cand:
            continue
        path = Path(cand).expanduser()
        if path.is_file() and os.access(path, os.X_OK):
            resolved = str(path.resolve())
            return f"{resolved} {rest}".strip()
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
    return problems
