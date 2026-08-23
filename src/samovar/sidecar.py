"""Optional sidecar conda envs for fragile tools (NanoSim, ART).

NanoSim 3.x needs Python 3.8 and runtime modules that are not fully declared
by every Bioconda build. It must not share the SamovaR env. ART is optional for
CAMISIM Illumina when Nextflow conda is disabled or you want a pre-built prefix.

Install::

    ./install.sh NanoSim
    ./install.sh ART

Prefixes default to ``~/.cache/samovar/envs/<name>`` and are recorded in
``tool_envs`` / ``nanosim_path`` / ``art_path``.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import List, Optional, Sequence

from samovar.paths import (
    conda_prefix_for_executable,
    discover_art,
    discover_nanosim,
    load_config,
    sidecar_envs_dir,
    write_config,
)

SIDECARS = {
    "nanosim": {
        "python": "3.8",
        "packages": [
            "python=3.8",
            "nanosim=3.2",
            "htseq",
            "gffutils",
            "numpy=1.23.5",
            "scikit-learn=0.23.2",
        ],
        "runtime_modules": (
            "HTSeq",
            "pysam",
            "sklearn.neighbors._dist_metrics",
        ),
        "runtime_versions": {"scikit-learn": "0.23.2"},
        "binary": "simulator.py",
        "config_path_key": "nanosim_path",
        "tool_keys": ("nanosim", "nanosim3", "simulator.py"),
        "hint": "ONT CAMISIM mode (./install.sh NanoSim). Do not install into the SamovaR env.",
    },
    "art": {
        "python": None,
        "packages": ["art", "samtools"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "art_illumina",
        "config_path_key": "art_path",
        "tool_keys": ("art", "art_illumina"),
        "hint": "Illumina CAMISIM mode (./install.sh ART).",
    },
}


def conda_executable() -> Optional[str]:
    env = os.environ.get("SAMOVAR_CONDA", "").strip()
    for name in (env, "mamba", "micromamba", "conda"):
        if not name:
            continue
        path = Path(name).expanduser()
        if path.is_file():
            return str(path.resolve())
        found = shutil.which(name)
        if found:
            return found
    return None


def sidecar_prefix(name: str) -> Path:
    return sidecar_envs_dir() / name


def env_has_binary(prefix: Path, binary: str) -> bool:
    cand = prefix / "bin" / binary
    try:
        return cand.is_file() and os.access(cand, os.X_OK)
    except OSError:
        return False


def env_has_runtime(prefix: Path, modules: Sequence[str]) -> bool:
    """Return true when the sidecar Python can import every required module."""
    if not modules:
        return True
    python = prefix / "bin" / "python"
    if not python.is_file():
        return False
    code = "; ".join(f"import {module}" for module in modules)
    try:
        result = subprocess.run(
            [str(python), "-c", code],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


def _purelib_path(prefix: Path) -> Optional[Path]:
    python = prefix / "bin" / "python"
    if not python.is_file():
        return None
    try:
        result = subprocess.run(
            [
                str(python),
                "-c",
                "import sysconfig; print(sysconfig.get_paths()['purelib'])",
            ],
            capture_output=True,
            text=True,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return None
    if result.returncode != 0 or not result.stdout.strip():
        return None
    return Path(result.stdout.strip().splitlines()[-1])


def remove_nanosim_compat_shim(prefix: Path) -> None:
    """Remove the temporary sklearn-1.x alias used by older local installs."""
    site_packages = _purelib_path(prefix)
    if site_packages is None:
        return
    try:
        (site_packages / "samovar_nanosim_sklearn_compat.pth").unlink(
            missing_ok=True
        )
    except OSError:
        return


def ensure_sidecar_compat(name: str, prefix: Path) -> None:
    if name == "nanosim":
        remove_nanosim_compat_shim(prefix)


def env_has_versions(prefix: Path, versions: dict) -> bool:
    if not versions:
        return True
    python = prefix / "bin" / "python"
    if not python.is_file():
        return False
    checks = " and ".join(
        f"version({package!r}) == {expected!r}"
        for package, expected in versions.items()
    )
    code = (
        "from importlib.metadata import version; "
        f"raise SystemExit(0 if ({checks}) else 1)"
    )
    try:
        result = subprocess.run(
            [str(python), "-c", code],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
            timeout=30,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return result.returncode == 0


def sidecar_is_healthy(name: str, prefix: Path) -> bool:
    spec = SIDECARS[name]
    return env_has_binary(prefix, str(spec["binary"])) and env_has_runtime(
        prefix, tuple(spec.get("runtime_modules") or ())
    ) and env_has_versions(prefix, dict(spec.get("runtime_versions") or {}))


def _conda_packages_cmd(
    conda: str,
    action: str,
    prefix: Path,
    packages: Sequence[str],
) -> List[str]:
    cmd = [conda, action, "-y", "-p", str(prefix)]
    cmd.extend(["-c", "conda-forge", "-c", "bioconda"])
    cmd.extend(packages)
    return cmd


def create_sidecar_env(
    name: str,
    *,
    packages: Optional[Sequence[str]] = None,
    python: Optional[str] = None,
) -> Path:
    spec = SIDECARS.get(name)
    if spec is None:
        raise ValueError(f"Unknown sidecar env {name!r}. Use: {', '.join(SIDECARS)}")
    prefix = sidecar_prefix(name)
    binary = str(spec["binary"])
    ensure_sidecar_compat(name, prefix)
    if sidecar_is_healthy(name, prefix):
        return prefix
    conda = conda_executable()
    if not conda:
        raise RuntimeError(
            "conda/mamba/micromamba is required to create a sidecar env. "
            f"Install {name} later with: ./install.sh {name.capitalize() if name != 'art' else 'ART'}"
        )
    pkgs = list(packages or spec["packages"])
    py = python if python is not None else spec.get("python")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    if py and not any(p.startswith("python") for p in pkgs):
        pkgs.insert(0, f"python={py}")
    action = "install" if prefix.is_dir() else "create"
    cmd = _conda_packages_cmd(conda, action, prefix, pkgs)
    verb = "Repairing" if action == "install" else "Creating"
    print(f"{verb} sidecar env:", " ".join(cmd))
    subprocess.check_call(cmd)
    ensure_sidecar_compat(name, prefix)
    if not sidecar_is_healthy(name, prefix):
        raise RuntimeError(
            f"Sidecar env {prefix} is incomplete after installation "
            f"({binary} or required runtime modules are missing)."
        )
    return prefix


def record_sidecar(name: str, prefix: Path, binary_path: str) -> None:
    spec = SIDECARS[name]
    cfg = load_config()
    cfg[str(spec["config_path_key"])] = binary_path
    tools = dict(cfg.get("tools") or {})
    for key in spec["tool_keys"]:
        tools[str(key)] = binary_path
    cfg["tools"] = tools
    envs = dict(cfg.get("tool_envs") or {})
    for key in spec["tool_keys"]:
        envs[str(key)] = str(prefix)
    cfg["tool_envs"] = envs
    write_config(cfg)


def install_named(name: str) -> Path:
    key = str(name or "").strip().lower()
    if key in {"art_illumina"}:
        key = "art"
    spec = SIDECARS.get(key)
    if spec is None:
        raise ValueError(f"Unknown sidecar {name!r}")
    existing = discover_nanosim() if key == "nanosim" else discover_art()
    if existing:
        prefix = Path(conda_prefix_for_executable(existing) or sidecar_prefix(key))
        ensure_sidecar_compat(key, prefix)
        if sidecar_is_healthy(key, prefix):
            record_sidecar(key, prefix, existing)
            print(f"{key} already present:", existing)
            print("tool_envs." + key + ":", prefix)
            return prefix
        print(f"{key} binary exists but its runtime is incomplete; repairing {prefix}")
    prefix = create_sidecar_env(key)
    binary = prefix / "bin" / str(spec["binary"])
    record_sidecar(key, prefix, str(binary.resolve()))
    print(f"{key}:", binary)
    print("tool_envs." + key + ":", prefix)
    return prefix


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Install a SamovaR sidecar conda env")
    parser.add_argument("name", choices=sorted(SIDECARS) + ["art_illumina"])
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    try:
        install_named(args.name)
    except (RuntimeError, subprocess.CalledProcessError) as exc:
        print("ERROR:", exc, file=sys.stderr)
        spec = SIDECARS.get("art" if args.name == "art_illumina" else args.name)
        if spec:
            print(spec["hint"], file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
