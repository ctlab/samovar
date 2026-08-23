"""Install-layout helpers: repo root, XDG config, and executable lookup.

Config is read from:

1. ``$SAMOVAR_CONFIG`` — exclusive override (even an empty ``{}``). Used by
   tests and custom installs so HPC ``~/.config`` / ``build/config.json``
   paths are not inherited.
2. Otherwise merge ``$XDG_CONFIG_HOME/samovar/config.json`` (default
   ``~/.config/samovar/config.json``) over ``<repo>/build/config.json``.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

PACKAGE_VERSION = "0.10.16"

KNOWN_TOOLS = (
    "kraken2",
    "kaiju",
    "kraken",
    "krakenuniq",
    "metaphlan",
    "metaphlan4",
    "centrifuge",
    "iss",
    "nextflow",
    "opal.py",
    "opal",
    "multiqc",
    "camisim",
    "simulator.py",
    "nanosim",
    "art_illumina",
    "wgsim",
    "samtools",
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
    if root:
        root_path = Path(root)
        try:
            usable = root_path.is_dir()
        except OSError:
            usable = False
        if usable:
            try:
                return root_path.resolve()
            except OSError:
                return root_path
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


def user_cache_dir() -> Path:
    xdg = os.environ.get("XDG_CACHE_HOME", "").strip()
    if xdg:
        return Path(xdg) / "samovar"
    return Path.home() / ".cache" / "samovar"


def genome_download_dir() -> Path:
    """Directory where SamovaR stores NCBI-downloaded (full) assemblies."""
    env = os.environ.get("SAMOVAR_GENOMES", "").strip()
    if env:
        return Path(env).expanduser()
    cfg = load_config()
    val = str(cfg.get("genomes") or "").strip()
    if val:
        return Path(val).expanduser()
    return user_cache_dir() / "genomes"


def processed_genomes_dir() -> Path:
    """Directory for ``{taxid}-processed.fasta.gz`` files (defaults to genomes cache)."""
    env = os.environ.get("SAMOVAR_PROCESSED_GENOMES", "").strip()
    if env:
        return Path(env).expanduser()
    cfg = load_config()
    val = str(cfg.get("processed_genomes") or "").strip()
    if val:
        return Path(val).expanduser()
    return genome_download_dir()


def update_config(updates: Dict[str, Any], also_repo_build: bool = True) -> Dict[str, Any]:
    """Merge ``updates`` into the install config and write it back."""
    cfg = load_config()
    changed = False
    for key, value in updates.items():
        if cfg.get(key) != value:
            cfg[key] = value
            changed = True
    if changed or not user_config_path().is_file():
        write_config(cfg, also_repo_build=also_repo_build)
    return cfg


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
    override = os.environ.get("SAMOVAR_CONFIG", "").strip()
    if override:
        # Exclusive: tests and custom installs must not inherit ~/.config or
        # build/config.json (those can point at unreadable HPC paths).
        data = _read_json(Path(override).expanduser())
        return dict(data) if isinstance(data, dict) else {}

    candidates = [
        user_config_dir() / "config.json",
        repo_root() / "build" / "config.json",
    ]
    merged: Dict[str, Any] = {}
    # Later files overlay earlier so repo build/ can fill gaps, user config wins.
    for path in reversed(list(dict.fromkeys(candidates))):
        data = _read_json(path)
        if data is not None:
            merged.update(data)
    return merged


def python_path() -> str:
    cfg = load_config()
    configured = (cfg.get("python_path") or "").strip()
    if configured:
        try:
            present = Path(configured).is_file()
        except OSError:
            present = False
        if present:
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
    try:
        is_file = path.is_file()
        is_dir = path.is_dir()
        bindir_ok = (path / "bin").is_dir()
    except OSError:
        is_file = False
        is_dir = False
        bindir_ok = False
    if is_file:
        return str(path.resolve().parent)
    bindir = path / "bin"
    if is_dir:
        if path.name == "bin":
            return str(path.resolve())
        if bindir_ok:
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
    add(cfg.get("opal_path"))
    add(cfg.get("multiqc_path"))
    add(cfg.get("nextflow_path"))
    add(cfg.get("nanosim_path"))
    add(cfg.get("art_path"))
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
    token_path = Path(token).expanduser()
    if (
        (token_path.is_absolute() or "/" in token)
        and token_path.is_file()
        and os.access(token_path, os.X_OK)
    ):
        return f"{str(token_path.resolve())} {rest}".strip()

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
    dest = user_config_path()
    dest.parent.mkdir(parents=True, exist_ok=True)
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
    opal = discover_opal()
    if opal:
        found.setdefault("opal.py", opal)
    multiqc = discover_multiqc()
    if multiqc:
        found.setdefault("multiqc", multiqc)
    try:
        from samovar.camisim import discover_camisim as _discover_camisim

        camisim = _discover_camisim()
        if camisim:
            found.setdefault("camisim", camisim)
    except Exception:
        pass
    nxt = shutil.which("nextflow")
    if nxt:
        found.setdefault("nextflow", nxt)
    nano = discover_nanosim()
    if nano:
        found.setdefault("simulator.py", nano)
        found.setdefault("nanosim", nano)
    art = discover_art()
    if art:
        found.setdefault("art_illumina", art)
    wgs = discover_wgsim()
    if wgs:
        found.setdefault("wgsim", wgs)
    return found


def discover_opal() -> Optional[str]:
    """Locate ``opal.py`` from env, config, PATH, or the install Python's scripts dir.

    ``cami-opal`` ships ``opal.py`` as a setuptools script (often under the
    interpreter's ``bin/`` or ``~/.local/bin``, which may not be on PATH).
    """
    cfg = load_config()
    candidates: List[str] = []

    def add(raw: Optional[str]) -> None:
        text = str(raw or "").strip()
        if text:
            candidates.append(text.split()[0])

    add(os.environ.get("SAMOVAR_OPAL_PATH") or os.environ.get("SAMOVAR_OPAL_BIN"))
    add(str(cfg.get("opal_path") or ""))
    tools = cfg.get("tools") if isinstance(cfg.get("tools"), dict) else {}
    add(str(tools.get("opal.py") or tools.get("opal") or ""))
    add(shutil.which("opal.py"))
    add(shutil.which("opal"))

    py = (cfg.get("python_path") or "").strip() or sys.executable
    py_path = Path(py).expanduser()
    if py_path.is_file():
        add(str(py_path.resolve().parent / "opal.py"))
        add(str(py_path.resolve().parent / "opal"))
    try:
        import sysconfig

        scripts = sysconfig.get_path("scripts")
        if scripts:
            add(str(Path(scripts) / "opal.py"))
            add(str(Path(scripts) / "opal"))
    except Exception:
        pass
    add(str(Path.home() / ".local" / "bin" / "opal.py"))
    add(str(Path.home() / ".local" / "bin" / "opal"))

    seen = set()
    for raw in candidates:
        if raw in seen:
            continue
        seen.add(raw)
        path = Path(raw).expanduser()
        if path.is_file():
            return str(path.resolve())
    return None


def discover_multiqc() -> Optional[str]:
    """Locate the MultiQC CLI (optional, same pattern as OPAL)."""
    cfg = load_config()
    candidates: List[str] = []

    def add(raw: Optional[str]) -> None:
        text = str(raw or "").strip()
        if text:
            candidates.append(text.split()[0])

    add(os.environ.get("SAMOVAR_MULTIQC_PATH") or os.environ.get("SAMOVAR_MULTIQC_BIN"))
    add(str(cfg.get("multiqc_path") or ""))
    tools = cfg.get("tools") if isinstance(cfg.get("tools"), dict) else {}
    add(str(tools.get("multiqc") or ""))
    add(shutil.which("multiqc"))

    py = (cfg.get("python_path") or "").strip() or sys.executable
    py_path = Path(py).expanduser()
    if py_path.is_file():
        add(str(py_path.resolve().parent / "multiqc"))
    try:
        import sysconfig

        scripts = sysconfig.get_path("scripts")
        if scripts:
            add(str(Path(scripts) / "multiqc"))
    except Exception:
        pass
    add(str(Path.home() / ".local" / "bin" / "multiqc"))

    seen = set()
    for raw in candidates:
        if raw in seen:
            continue
        seen.add(raw)
        path = Path(raw).expanduser()
        if path.is_file():
            return str(path.resolve())
    return None


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


def sidecar_envs_dir() -> Path:
    """Conda prefixes for fragile optional tools (NanoSim, ART, …)."""
    return user_cache_dir() / "envs"


def tool_env_prefix(name: str, cfg: Optional[Dict[str, Any]] = None) -> Optional[str]:
    """Absolute conda/module prefix from ``tool_envs.<name>`` if it exists."""
    cfg = dict(cfg or load_config())
    envs = cfg.get("tool_envs") if isinstance(cfg.get("tool_envs"), dict) else {}
    aliases = (name, name.replace("-", "_"), Path(name).stem)
    for key in aliases:
        raw = str(envs.get(key) or "").strip()
        if not raw:
            continue
        path = Path(raw).expanduser()
        try:
            if path.is_dir():
                return str(path.resolve())
        except OSError:
            continue
        return str(path)
    return None


def _first_existing_executable(candidates: List[str]) -> Optional[str]:
    seen = set()
    for raw in candidates:
        text = str(raw or "").strip()
        if not text or text in seen:
            continue
        seen.add(text)
        path = Path(text).expanduser()
        try:
            if path.is_file() and os.access(path, os.X_OK):
                return str(path.resolve())
        except OSError:
            continue
        found = shutil.which(text)
        if found:
            return found
    return None


def discover_nanosim() -> Optional[str]:
    """Locate NanoSim ``simulator.py`` (optional; keep it out of the main env)."""
    cfg = load_config()
    candidates: List[str] = []

    def add(raw: Optional[str]) -> None:
        text = str(raw or "").strip()
        if text:
            candidates.append(text.split()[0])

    add(os.environ.get("SAMOVAR_NANOSIM") or os.environ.get("SAMOVAR_NANOSIM_BIN"))
    add(str(cfg.get("nanosim_path") or ""))
    tools = cfg.get("tools") if isinstance(cfg.get("tools"), dict) else {}
    add(str(tools.get("simulator.py") or tools.get("nanosim") or tools.get("nanosim3") or ""))
    prefix = tool_env_prefix("nanosim", cfg) or tool_env_prefix("nanosim3", cfg)
    if prefix:
        add(str(Path(prefix) / "bin" / "simulator.py"))
        add(str(Path(prefix) / "bin" / "nanosim"))
    default_prefix = sidecar_envs_dir() / "nanosim"
    add(str(default_prefix / "bin" / "simulator.py"))
    add(shutil.which("simulator.py"))
    add(shutil.which("nanosim"))
    return _first_existing_executable(candidates)


def discover_art() -> Optional[str]:
    """Locate ``art_illumina`` (optional ART Illumina simulator)."""
    cfg = load_config()
    candidates: List[str] = []

    def add(raw: Optional[str]) -> None:
        text = str(raw or "").strip()
        if text:
            candidates.append(text.split()[0])

    add(os.environ.get("SAMOVAR_ART") or os.environ.get("SAMOVAR_ART_ILLUMINA"))
    add(str(cfg.get("art_path") or ""))
    tools = cfg.get("tools") if isinstance(cfg.get("tools"), dict) else {}
    add(str(tools.get("art_illumina") or tools.get("art") or ""))
    prefix = tool_env_prefix("art", cfg) or tool_env_prefix("art_illumina", cfg)
    if prefix:
        add(str(Path(prefix) / "bin" / "art_illumina"))
    add(str(sidecar_envs_dir() / "art" / "bin" / "art_illumina"))
    add(shutil.which("art_illumina"))
    return _first_existing_executable(candidates)


def discover_wgsim() -> Optional[str]:
    cfg = load_config()
    candidates: List[str] = []

    def add(raw: Optional[str]) -> None:
        text = str(raw or "").strip()
        if text:
            candidates.append(text.split()[0])

    add(os.environ.get("SAMOVAR_WGSIM"))
    add(str(cfg.get("wgsim_path") or ""))
    tools = cfg.get("tools") if isinstance(cfg.get("tools"), dict) else {}
    add(str(tools.get("wgsim") or ""))
    prefix = tool_env_prefix("wgsim", cfg)
    if prefix:
        add(str(Path(prefix) / "bin" / "wgsim"))
    add(shutil.which("wgsim"))
    return _first_existing_executable(candidates)


def conda_prefix_for_executable(exe: Optional[str]) -> Optional[str]:
    """If ``exe`` lives in ``<prefix>/bin``, return that conda-style prefix."""
    if not exe:
        return None
    path = Path(exe).expanduser()
    try:
        if not path.is_file():
            return None
        parent = path.resolve().parent
    except OSError:
        return None
    if parent.name == "bin" and parent.parent.is_dir():
        return str(parent.parent)
    return None


def _existing_tool(name: str) -> Optional[str]:
    resolved = resolve_executable(name, tool_key=name)
    token = str(resolved or "").split()[0]
    if not token:
        return None
    path = Path(token)
    try:
        if path.is_file() and os.access(path, os.X_OK):
            return str(path.resolve())
    except OSError:
        pass
    return shutil.which(name)


def install_status_rows() -> List[Dict[str, Any]]:
    """Required and optional tools: path, role, how to install, config key."""
    cfg = load_config()
    iss = iss_executable()
    iss_ok = bool(iss and (Path(iss).is_file() or shutil.which(str(iss).split()[0])))
    combiner = repo_root() / "bin" / "samovar_combine_annotations"
    viz = ""
    try:
        from samovar.viz_annotation import require_viz_backend

        viz = require_viz_backend()
    except Exception:
        viz = ""

    def row(
        name: str,
        *,
        required: bool,
        role: str,
        path: Optional[str],
        config_key: str,
        install: str,
    ) -> Dict[str, Any]:
        found = bool(path)
        return {
            "name": name,
            "required": required,
            "role": role,
            "path": path or "",
            "found": found,
            "config_key": config_key,
            "install": install,
        }

    py = python_path()
    rows = [
        row(
            "python",
            required=True,
            role="SamovaR runtime",
            path=py if Path(py).is_file() else shutil.which("python3"),
            config_key="python_path",
            install="conda env create -f environment.yml",
        ),
        row(
            "iss",
            required=True,
            role="InSilicoSeq read simulation (default generate)",
            path=iss if iss_ok else None,
            config_key="iss_path",
            install="pip install 'insilicoseq>2.0.0'  (or conda: insilicoseq)",
        ),
        row(
            "snakemake",
            required=True,
            role="Annotator workflow runner",
            path=shutil.which("snakemake"),
            config_key="tools.snakemake",
            install="pip install 'snakemake>=7.0'  (or conda: snakemake-minimal)",
        ),
        row(
            "C++ combiner",
            required=True,
            role="Merge annotator reports",
            path=str(combiner) if combiner.is_file() else cxx_compiler(),
            config_key="",
            install="g++ / make -C src/cpp  (install.sh builds this)",
        ),
        row(
            "viz backend",
            required=True,
            role="Pipeline plots (cnsplots or altair)",
            path=viz or None,
            config_key="",
            install="pip install 'cnsplots>=0.6.0' altair",
        ),
        row(
            "kraken2",
            required=False,
            role="k-mer LCA annotator",
            path=_existing_tool("kraken2"),
            config_key="tools.kraken2 / tool_envs.kraken2",
            install="conda create -c bioconda -n kraken2 kraken2; set tool_envs.kraken2",
        ),
        row(
            "kaiju",
            required=False,
            role="protein-level annotator",
            path=_existing_tool("kaiju"),
            config_key="tools.kaiju / tool_envs.kaiju",
            install="conda create -c bioconda -n kaiju kaiju; set tool_envs.kaiju",
        ),
        row(
            "metaphlan",
            required=False,
            role="marker-gene profiler",
            path=shutil.which("metaphlan") or shutil.which("metaphlan4"),
            config_key="tools.metaphlan / tool_envs.metaphlan",
            install="conda install -c bioconda metaphlan",
        ),
        row(
            "centrifuge",
            required=False,
            role="optional annotator",
            path=shutil.which("centrifuge"),
            config_key="tools.centrifuge / tool_envs.centrifuge",
            install="conda install -c bioconda centrifuge",
        ),
        row(
            "OPAL",
            required=False,
            role="CAMI HTML profile assessment (metrics always computed in Python)",
            path=discover_opal(),
            config_key="opal_path",
            install="./install.sh OPAL",
        ),
        row(
            "MultiQC",
            required=False,
            role="HTML report of SamovaR plots",
            path=discover_multiqc(),
            config_key="multiqc_path",
            install="./install.sh MultiQC",
        ),
        row(
            "CAMISIM",
            required=False,
            role="Optional generate backend (community design / ART / NanoSim / wgsim)",
            path=None,
            config_key="camisim_path",
            install="./install.sh CAMISIM",
        ),
        row(
            "Nextflow",
            required=False,
            role="CAMISIM 2.0 runner (table mode does not need it)",
            path=shutil.which("nextflow") or str(cfg.get("nextflow_path") or "") or None,
            config_key="nextflow_path / tool_envs.nextflow",
            install="conda install -c bioconda nextflow   or https://nextflow.io",
        ),
        row(
            "NanoSim",
            required=False,
            role="ONT reads for CAMISIM --camisim-mode ont|hybrid (separate env)",
            path=discover_nanosim(),
            config_key="nanosim_path / tool_envs.nanosim",
            install="./install.sh NanoSim",
        ),
        row(
            "ART",
            required=False,
            role="Illumina reads for CAMISIM --camisim-mode illumina|hybrid",
            path=discover_art(),
            config_key="art_path / tool_envs.art",
            install="./install.sh ART",
        ),
        row(
            "wgsim",
            required=False,
            role="Simple paired reads for CAMISIM --camisim-mode wgsim",
            path=discover_wgsim(),
            config_key="wgsim_path / tool_envs.wgsim",
            install="conda install -c bioconda wgsim   (or ./install.sh ART which can share samtools)",
        ),
        row(
            "R",
            required=False,
            role="Optional samovaR regenerator (GitHub r-package branch)",
            path=shutil.which("R") or str(cfg.get("r_path") or "") or None,
            config_key="r_path",
            install="./install.sh R-package",
        ),
    ]
    try:
        from samovar.camisim import discover_camisim as _dc

        cam = _dc()
    except Exception:
        cam = str(cfg.get("camisim_path") or "") or None
    for item in rows:
        if item["name"] == "CAMISIM":
            item["path"] = cam or ""
            item["found"] = bool(cam)
        if item["name"] == "Nextflow":
            nxt = str(cfg.get("nextflow_path") or "").strip() or shutil.which("nextflow")
            token = str(nxt or "").split()[0]
            ok = bool(token and (Path(token).is_file() or shutil.which("nextflow")))
            item["path"] = (str(Path(token).resolve()) if token and Path(token).is_file() else shutil.which("nextflow")) or ""
            item["found"] = ok and bool(item["path"])
    return rows


def format_install_status() -> str:
    """Human-readable required/optional tool table for ``install.sh`` and ``samovar tools --status``."""
    rows = install_status_rows()
    lines = [
        "SamovaR tool status",
        "===================",
        "Required:",
    ]
    for item in rows:
        if not item["required"]:
            continue
        mark = "ok" if item["found"] else "MISSING"
        loc = item["path"] or item["install"]
        lines.append(f"  [{mark:7}] {item['name']:<14} {loc}")
    lines.append("Optional:")
    for item in rows:
        if item["required"]:
            continue
        mark = "ok" if item["found"] else "—"
        loc = item["path"] if item["found"] else item["install"]
        lines.append(f"  [{mark:7}] {item['name']:<14} {loc}")
        if not item["found"]:
            lines.append(f"           {item['role']}")
            if item["config_key"]:
                lines.append(f"           config: {item['config_key']}")
    missing = [i["name"] for i in rows if i["required"] and not i["found"]]
    if missing:
        lines.append("")
        lines.append("Required tools missing: " + ", ".join(missing))
    else:
        lines.append("")
        lines.append("All required programs are present.")
    return "\n".join(lines)
