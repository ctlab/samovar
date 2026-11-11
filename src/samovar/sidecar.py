"""Optional sidecar conda envs for fragile tools (NanoSim, ART).

NanoSim 3.x needs Python 3.8 and runtime modules that are not fully declared
by every Bioconda build. It must not share the SamovaR env. ART is optional for
CAMISIM Illumina when Nextflow conda is disabled or you want a pre-built prefix.

Install::

    ./install.sh NanoSim
    ./install.sh ART

Prefixes default to ``$SAMOVAR_ROOT/.cache/samovar/envs/<name>`` (override with
``SAMOVAR_ENVS``) and are recorded in ``tool_envs`` / ``nanosim_path`` /
``art_path``. Never under ``~/.cache`` by default (home quota).
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
from samovar.main_config import lookup_tool_record, set_tool

SIDECARS = {
    "phanta": {
        "python": None, "packages": ["kraken2>=2.1.6"], "runtime_modules": (), "runtime_versions": {},
        "binary": "kraken2", "tool_keys": ("phanta",), "group": "profiler",
        "required_help": "--report-minimizer-data",
        "dependencies": ("kraken2", "bracken"),
        "source": {
            "url": "https://codeload.github.com/bhattlab/phanta/tar.gz/ecbfcf9fbf874072ec6069913291991cbec6bab3",
            "sha256": "e1a2efe568b281f23e6f81eb3075481437044acbcecac1b99088c7267cd3c115",
            "revision": "ecbfcf9fbf874072ec6069913291991cbec6bab3",
            "directory": "phanta-ecbfcf9fbf874072ec6069913291991cbec6bab3",
            "executable": "run_phanta.py",
        },
        "hint": "Phanta native sample profiler (./install.sh Phanta); database installed separately.",
    },
    "motus": {
        "python": "3.8", "packages": ["python=3.8", "python_abi=3.8=*_cp38", "motus=3.1.0", "bwa", "samtools"], "runtime_modules": (),
        "runtime_versions": {"motu": "3.1.0"}, "binary": "motus", "tool_keys": ("motus",),
        "companions": ("bwa", "samtools"), "companion_group": "aligner", "group": "profiler",
        "hint": "mOTUs3 sample-level profiler (./install.sh mOTUs3); reference DB installed separately.",
    },
    "kmcp": {
        "python": None, "packages": ["kmcp>=0.9.5"], "runtime_modules": (),
        "runtime_versions": {}, "binary": "kmcp", "tool_keys": ("kmcp",),
        "group": "annotator", "hint": "Native KMCP read annotator (./install.sh KMCP).",
    },
    "bracken": {
        "python": None, "packages": ["bracken"], "runtime_modules": (),
        "runtime_versions": {}, "binary": "bracken", "tool_keys": ("bracken",),
        "companions": ("bracken-build",), "group": "profiler",
        "hint": "Kraken2 + Bracken sample-level profiler (./install.sh Bracken).",
    },
    "sylph": {
        "python": None, "packages": ["sylph>=0.6"], "runtime_modules": (),
        "runtime_versions": {}, "binary": "sylph", "tool_keys": ("sylph",),
        "group": "profiler", "hint": "Optional sample-level profiler (./install.sh Sylph).",
    },
    "kraken2": {
        "python": None, "packages": ["kraken2"], "runtime_modules": (),
        "runtime_versions": {}, "binary": "kraken2", "tool_keys": ("kraken2",),
        "companions": ("kraken2-build", "kraken2-inspect"), "group": "annotator",
        "hint": "Core ensemble annotator (./install.sh Kraken2; included in full).",
    },
    "kaiju": {
        "python": None, "packages": ["kaiju"], "runtime_modules": (),
        "runtime_versions": {}, "binary": "kaiju", "tool_keys": ("kaiju",),
        "companions": ("kaiju-mkbwt", "kaiju-mkfmi", "kaiju-multi"), "group": "annotator",
        "hint": "Core ensemble annotator (./install.sh Kaiju; included in full).",
    },
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
    "megahit": {
        "python": None,
        "packages": ["megahit"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "megahit",
        "config_path_key": "megahit_path",
        "tool_keys": ("megahit",),
        "group": "assembler",
        "hint": "Optional MegaHIT assembler (./install.sh MegaHIT).",
    },
    "flye": {
        "python": None,
        "packages": ["flye"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "flye",
        "config_path_key": "flye_path",
        "tool_keys": ("flye",),
        "group": "assembler",
        "hint": "Optional Flye / metaFlye assembler (./install.sh metaFlye).",
    },
    "spades": {
        "python": None,
        "packages": ["spades"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "spades.py",
        "config_path_key": "spades_path",
        "tool_keys": ("spades",),
        "group": "assembler",
        "hint": "Optional SPAdes / metaSPAdes assembler (./install.sh SPAdes; included in full).",
    },
    "metax": {
        "python": None,
        "packages": ["metax=0.9.23=beta_0", "ma=1.1.4"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "metax",
        "extra_channels": ("zldeng",),
        "config_path_key": "metax_path",
        "tool_keys": ("metax",),
        "group": "annotator",
        "hint": "Optional Metax 0.9.23 profiler (./install.sh Metax). Pinned below 2.x so prebuilt RefSeq indexes (metax_db.json + ma) still load. Needs maCMD from bioconda ma=1.1.4.",
    },
    "prodigal": {
        "python": None,
        "packages": ["prodigal"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "prodigal",
        "config_path_key": "prodigal_path",
        "tool_keys": ("prodigal",),
        "group": "gene_caller",
        "hint": "Optional Prodigal gene caller (./install.sh Prodigal).",
    },
    "minimap2": {
        "python": None,
        "packages": ["minimap2", "samtools"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "minimap2",
        "config_path_key": "minimap2_path",
        "tool_keys": ("minimap2",),
        "group": "aligner",
        "hint": "Optional minimap2 MAG aligner (./install.sh minimap2).",
    },
    "coverm": {
        "python": None,
        "packages": ["coverm"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "coverm",
        "config_path_key": "coverm_path",
        "tool_keys": ("coverm",),
        "group": "mag_quantifier",
        "hint": "Optional CoverM MAG quantifier (./install.sh CoverM).",
    },
    "dastool": {
        "python": None,
        "packages": ["das_tool"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "DAS_Tool",
        "config_path_key": "dastool_path",
        "tool_keys": ("dastool", "DAS_Tool"),
        "group": "binner_combine",
        "hint": "Optional DAS Tool MAG combine (./install.sh DAS_Tool).",
    },
    "anvio": {
        "python": "3.10",
        "packages": [
            "python=3.10",
            "anvio-minimal=8",
            "concoct",
            "prodigal",
            "samtools",
            "minimap2",
            "setuptools",
        ],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "anvi-gen-contigs-database",
        "config_path_key": "anvio_path",
        "tool_keys": ("anvio", "anvi-gen-contigs-database"),
        "group": "binner",
        "hint": "Optional anvi'o binner sidecar (./install.sh anvio).",
    },
    "checkm2": {
        "python": None,
        "packages": ["checkm2"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "checkm2",
        "config_path_key": "checkm2_path",
        "tool_keys": ("checkm2",),
        "group": "binner_qc",
        "hint": "Optional CheckM2 sidecar (./install.sh CheckM2).",
    },
    "gtdbtk": {
        "python": "3.10",
        "packages": ["python=3.10", "gtdbtk=2.7.2", "prodigal"],
        "runtime_modules": (),
        "runtime_versions": {},
        "binary": "gtdbtk",
        "config_path_key": "gtdbtk_path",
        "tool_keys": ("gtdbtk",),
        "group": "mag_taxonomy",
        "hint": "Optional GTDB-Tk sidecar (./install.sh GTDB-Tk). Set GTDBTK_DATA_PATH.",
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
    if name == "anvio":
        python = prefix / "bin" / "python"
        if python.is_file():
            subprocess.run(
                [str(python), "-m", "pip", "install", "setuptools<81"],
                check=False,
                timeout=180,
            )


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
    if spec.get('required_help'):
        binary = prefix / 'bin' / spec['binary']
        if not binary.is_file():
            return False
        result = subprocess.run([str(binary), '--help'], capture_output=True, text=True, timeout=30)
        if spec['required_help'] not in result.stdout + result.stderr:
            return False
    return all(env_has_binary(prefix, binary) for binary in (str(spec["binary"]), *spec.get("companions", ()))) and env_has_runtime(
        prefix, tuple(spec.get("runtime_modules") or ())
    ) and env_has_versions(prefix, dict(spec.get("runtime_versions") or {}))


def _conda_packages_cmd(
    conda: str,
    action: str,
    prefix: Path,
    packages: Sequence[str],
    extra_channels: Sequence[str] = (),
) -> List[str]:
    cmd = [conda, action, "-y", "-p", str(prefix)]
    cmd.extend(["-c", "conda-forge", "-c", "bioconda"])
    for channel in extra_channels:
        cmd.extend(["-c", str(channel)])
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
    cmd = _conda_packages_cmd(
        conda,
        action,
        prefix,
        pkgs,
        extra_channels=tuple(spec.get("extra_channels") or ()),
    )
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
    from samovar.main_config import LEGACY_TOOL_KEYS

    path_key = str(spec.get("config_path_key") or "")
    if path_key in LEGACY_TOOL_KEYS:
        cfg[path_key] = binary_path
    for key in spec["tool_keys"]:
        group = str(spec.get("group") or (
            "metagenome_generator" if name == "nanosim" else "reads_generator"
        ))
        set_tool(
            cfg,
            str(key),
            env="conda",
            workflow=str(name),
            path=str(prefix),
            group=group,
        )
        # Keep the runnable binary on the primary name
        if str(key) in {spec["binary"], name}:
            set_tool(
                cfg,
                str(key),
                env="conda",
                workflow=str(name),
                path=binary_path,
                group=group,
            )
    for companion in spec.get("companions", ()):
        if lookup_tool_record(cfg, companion):
            continue  # Keep independently configured companion programs.
        set_tool(cfg, companion, env="conda", workflow=name,
                 path=str((prefix / "bin" / companion).resolve()),
                 group=str(spec.get("companion_group") or spec.get("group") or "annotator"))
    write_config(cfg)


def _install_source_tool(name: str) -> Path:
    """Install a pinned source workflow through the same optional installer."""
    import hashlib
    import tarfile
    import tempfile
    import urllib.request
    from samovar.paths import _existing_tool

    spec = SIDECARS[name]
    source = spec['source']
    for dependency in spec.get('dependencies', ()):
        if not _existing_tool(dependency):
            install_named(dependency)
    dest = sidecar_envs_dir().parent / 'sources' / name / source['revision']
    executable = dest / source['executable']
    if not executable.is_file():
        if os.environ.get('SAMOVAR_OFFLINE', '0') != '0':
            raise RuntimeError(f'{name} source is unavailable in offline mode: {dest}')
        dest.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=dest.parent, prefix='.install-') as tmp:
            work = Path(tmp)
            archive = work / 'source.tar.gz'
            with urllib.request.urlopen(source['url'], timeout=120) as response, archive.open('wb') as handle:
                shutil.copyfileobj(response, handle)
            if hashlib.sha256(archive.read_bytes()).hexdigest() != source['sha256']:
                raise RuntimeError(f'{name} source checksum mismatch')
            with tarfile.open(archive) as contents:
                kwargs = {'filter': 'data'} if hasattr(tarfile, 'data_filter') else {}
                contents.extractall(work / 'unpacked', **kwargs)
            (work / 'unpacked' / source['directory']).replace(dest)
    executable.chmod(executable.stat().st_mode | 0o111)
    cfg = load_config()
    kraken = _existing_tool('phanta-kraken2') or _existing_tool('kraken2')
    probe = subprocess.run([kraken, '--help'], capture_output=True, text=True, timeout=30)
    if spec['required_help'] not in probe.stdout + probe.stderr:
        prefix = create_sidecar_env(name)
        set_tool(cfg, 'phanta-kraken2', path=str(prefix / 'bin/kraken2'), group='annotator')
    set_tool(cfg, name, path=str(executable), group=spec['group'], version=source['revision'])
    write_config(cfg)
    print(f'{name}: {executable} (source revision {source["revision"]})')
    return dest


def install_named(name: str) -> Path:
    key = str(name or "").strip().lower()
    aliases = {
        "motus3": "motus",
        "art_illumina": "art",
        "das_tool": "dastool",
        "das-tool": "dastool",
        "gtdb-tk": "gtdbtk",
        "gtdb_tk": "gtdbtk",
        "mega-hit": "megahit",
        "metaflye": "flye",
        "metaspades": "spades",
        "meta-spades": "spades",
        "metaspades.py": "spades",
        "metax": "metax",
        "meta-flye": "flye",
        "metflye": "flye",
        "cover-m": "coverm",
    }
    key = aliases.get(key, key)
    spec = SIDECARS.get(key)
    if spec is None:
        raise ValueError(f"Unknown sidecar {name!r}")
    if spec.get('source'):
        return _install_source_tool(key)
    existing = None
    if key == "nanosim":
        existing = discover_nanosim()
    elif key == "art":
        existing = discover_art()
    else:
        if key in {"kraken2", "kaiju"}:
            from samovar.paths import _existing_tool
            existing = _existing_tool(key)
        else:
            existing = shutil.which(str(spec["binary"]))
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
    parser.add_argument("name", choices=sorted(SIDECARS) + ["art_illumina", "DAS_Tool", "GTDB-Tk", "MegaHIT", "CoverM", "CheckM2", "Prodigal"])
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
