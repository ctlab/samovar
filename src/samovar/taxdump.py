"""Shared NCBI taxdump location (nodes.dmp / names.dmp).

Configured as ``genomes.taxdump`` (env ``SAMOVAR_TAXDUMP``). Default:
``{samovar_database}/taxdump``. Downloads land there once; Kaiju/Kraken DBs
get symlinks to the ``.dmp`` files instead of a fresh NCBI fetch.
"""

from __future__ import annotations

import logging
import os
import tarfile
from pathlib import Path
from typing import Any, Dict, Iterable, Optional, Union

logger = logging.getLogger(__name__)

PathLike = Union[str, os.PathLike]

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Do not link prot.accession2taxid or other bulky extras into annotator DBs.
TAXDUMP_DMP_FILES = (
    "nodes.dmp",
    "names.dmp",
    "merged.dmp",
    "delnodes.dmp",
    "division.dmp",
    "gencode.dmp",
)


def _as_path(raw: PathLike) -> Path:
    return Path(raw).expanduser()


def taxdump_dir(cfg: Optional[Dict[str, Any]] = None) -> Path:
    """Install-level NCBI taxdump directory."""
    env = os.environ.get("SAMOVAR_TAXDUMP", "").strip()
    if env:
        return _as_path(env)
    if cfg is None:
        try:
            from samovar.paths import load_config

            cfg = load_config()
        except Exception:
            cfg = {}
    from samovar.main_config import genomes_block
    from samovar.genome_index import samovar_database_dir

    block = genomes_block(cfg or {}) if cfg else {}
    raw = str((block.get("taxdump") if isinstance(block, dict) else "") or "").strip()
    if raw:
        return _as_path(raw)
    return samovar_database_dir(cfg) / "taxdump"


def find_dmp(name: str, root: Optional[PathLike] = None) -> Optional[Path]:
    """``nodes.dmp`` / ``names.dmp`` under a taxdump dir (or ``taxonomy/``)."""
    base = _as_path(root) if root is not None else taxdump_dir()
    for candidate in (base / name, base / "taxonomy" / name):
        try:
            if candidate.is_file():
                return candidate
        except OSError:
            continue
    return None


def nodes_dmp(cfg: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    env = os.environ.get("SAMOVAR_NODES_DMP", "").strip()
    if env:
        path = _as_path(env)
        if path.is_file():
            return path
    return find_dmp("nodes.dmp", taxdump_dir(cfg))


def names_dmp(cfg: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    return find_dmp("names.dmp", taxdump_dir(cfg))


def taxdump_tarball(cfg: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    root = taxdump_dir(cfg)
    for name in ("taxdump.tar.gz", "ncbi-taxonomy.tar.gz"):
        path = root / name
        try:
            if path.is_file():
                return path
        except OSError:
            continue
    return None


def _extract_dmp_from_tarball(archive: Path, dest: Path) -> None:
    dest.mkdir(parents=True, exist_ok=True)
    wanted = set(TAXDUMP_DMP_FILES)
    with tarfile.open(archive, "r:*") as tar:
        for info in tar.getmembers():
            name = Path(info.name.replace("\\", "/")).name
            if name not in wanted or not info.isfile():
                continue
            target = dest / name
            handle = tar.extractfile(info)
            if handle is None:
                continue
            target.write_bytes(handle.read())
            logger.info("Extracted %s -> %s", name, target)


def ensure_taxdump(cfg: Optional[Dict[str, Any]] = None) -> Path:
    """Return a directory that contains ``nodes.dmp``, downloading once if needed."""
    dest = taxdump_dir(cfg)
    if find_dmp("nodes.dmp", dest) is not None:
        return dest
    dest.mkdir(parents=True, exist_ok=True)
    archive = taxdump_tarball(cfg)
    if archive is None:
        url = os.environ.get("SAMOVAR_TAXDUMP_URL", "").strip() or TAXDUMP_URL
        archive = dest / "taxdump.tar.gz"
        logger.info("Downloading NCBI taxdump to %s", archive)
        import subprocess

        subprocess.run(["wget", "-O", str(archive), url], check=True)
    logger.info("Extracting taxdump .dmp files into %s", dest)
    _extract_dmp_from_tarball(archive, dest)
    if find_dmp("nodes.dmp", dest) is None:
        raise FileNotFoundError(f"nodes.dmp missing after taxdump extract in {dest}")
    return dest


def _replace_symlink(dest: Path, source: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    try:
        if dest.is_symlink() or dest.is_file():
            dest.unlink()
        elif dest.exists():
            return
    except OSError:
        return
    os.symlink(source, dest)


def link_taxdump_into(db_path: PathLike, source: Optional[PathLike] = None) -> Path:
    """Symlink ``nodes.dmp`` / ``names.dmp`` (and siblings) into ``db_path``."""
    src_root = _as_path(source) if source else ensure_taxdump()
    if src_root.is_file() and str(src_root).endswith(".gz"):
        src_root = ensure_taxdump()
    dest = _as_path(db_path)
    dest.mkdir(parents=True, exist_ok=True)
    linked = 0
    for name in TAXDUMP_DMP_FILES:
        found = find_dmp(name, src_root)
        if found is None:
            continue
        _replace_symlink(dest / name, found.resolve())
        linked += 1
    if linked == 0:
        raise FileNotFoundError(f"No taxdump .dmp files under {src_root}")
    logger.info("Linked %s taxdump file(s) into %s", linked, dest)
    return dest


def dmp_search_roots(cfg: Optional[Dict[str, Any]] = None) -> Iterable[str]:
    """Directories/files to scan for ``nodes.dmp`` (env + config)."""
    roots = []
    env_file = os.environ.get("SAMOVAR_NODES_DMP", "").strip()
    if env_file:
        roots.append(env_file)
    extra = os.environ.get("SAMOVAR_NODES_SEARCH", "")
    roots.extend(part.strip() for part in extra.split(":") if part.strip())
    try:
        dump = taxdump_dir(cfg)
        roots.append(str(dump))
        roots.append(str(dump / "taxonomy"))
    except Exception:
        pass
    return roots
