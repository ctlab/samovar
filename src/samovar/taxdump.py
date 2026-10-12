"""Shared NCBI taxdump location (nodes.dmp / names.dmp).

Configured as ``genomes.taxdump`` (env ``SAMOVAR_TAXDUMP``). Default:
``{samovar_database}/taxdump``. Downloads land there once; Kaiju/Kraken DBs
get symlinks to the ``.dmp`` files instead of a fresh NCBI fetch.
"""

from __future__ import annotations

import hashlib
import logging
import os
import pickle
import tarfile
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Tuple, Union

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


def taxdump_version_label(root: PathLike) -> str:
    """Version stamp for one dump: ``taxdump.version``, else ``nodes.dmp`` date.

    NCBI does not ship a version file. The ``nodes.dmp`` mtime is the date the
    dump was unpacked, which is enough to keep two dumps as separate catalog
    keys (``ncbi:YYYY-MM-DD``).
    """
    base = _as_path(root)
    for name in ("taxdump.version", "VERSION"):
        path = base / name
        try:
            if path.is_file():
                text = path.read_text(encoding="utf-8", errors="ignore").strip().split()
                if text:
                    return text[0]
        except OSError:
            continue
    nodes = find_dmp("nodes.dmp", base)
    if nodes is None:
        return ""
    try:
        return datetime.fromtimestamp(nodes.stat().st_mtime).strftime("%Y-%m-%d")
    except OSError:
        return ""


def _record_path(rec: Dict[str, Any]) -> Optional[Path]:
    raw = str(rec.get("path") or "").strip()
    if not raw:
        return None
    return _as_path(raw)


def configured_taxdumps(
    cfg: Optional[Dict[str, Any]] = None,
) -> List[Tuple[str, Path]]:
    """``(version, directory)`` from env, ``genomes.taxdump``, and the catalog.

    Version may be empty. Directories are not required to exist.
    """
    if cfg is None:
        try:
            from samovar.paths import load_config

            cfg = load_config()
        except Exception:
            cfg = {}
    found: List[Tuple[str, Path]] = []
    seen = set()

    def add(version: str, path: Optional[Path]) -> None:
        if path is None:
            return
        key = str(path)
        if key in seen:
            return
        seen.add(key)
        found.append((str(version or "").strip(), path))

    env = os.environ.get("SAMOVAR_TAXDUMP", "").strip()
    if env:
        add(os.environ.get("SAMOVAR_TAXDUMP_VERSION", "").strip(), _as_path(env))
    from samovar.main_config import genomes_block

    block = genomes_block(cfg or {}) if cfg else {}
    if isinstance(block, dict):
        raw = str(block.get("taxdump") or "").strip()
        if raw:
            add(str(block.get("taxdump_version") or "").strip(), _as_path(raw))
    try:
        from samovar.db_spec import iter_database_records, split_db_key

        grouped = iter_database_records(cfg or {}).get("taxdump") or {}
        for key, item in grouped.items():
            if not isinstance(item, dict):
                continue
            bare, key_ver = split_db_key(str(key))
            if bare in {"gtdb", "gtdb_taxdump"}:
                continue
            version = str(item.get("version") or item.get("_version") or key_ver or "").strip()
            add(version, _record_path(item))
    except Exception:
        pass
    return found


def usable_taxdump_dir(
    cfg: Optional[Dict[str, Any]] = None,
    version: str = "",
) -> Optional[Path]:
    """A configured dump that already contains ``nodes.dmp``. Never downloads.

    ``version`` (or ``SAMOVAR_TAXDUMP_VERSION``) selects one stamp. Without it,
    ``genomes.taxdump`` wins when that directory is populated, otherwise the
    newest dated catalog entry.
    """
    want = str(version or os.environ.get("SAMOVAR_TAXDUMP_VERSION", "") or "").strip()
    live: List[Tuple[str, Path]] = []
    for stamp, path in configured_taxdumps(cfg):
        if find_dmp("nodes.dmp", path) is None:
            continue
        live.append((stamp or taxdump_version_label(path), path))
    if not live:
        return None
    if want:
        for stamp, path in live:
            if stamp == want:
                return path
        return None
    env = os.environ.get("SAMOVAR_TAXDUMP", "").strip()
    if env:
        env_path = _as_path(env)
        for _stamp, path in live:
            if path == env_path:
                return path
    return live[0][1]


def metax_dmp_dir(cfg: Optional[Dict[str, Any]] = None, version: str = "") -> str:
    """Directory Metax should pass as ``--dmp-dir``. Empty if none is installed."""
    found = usable_taxdump_dir(cfg, version=version)
    if found is None:
        return ""
    nodes = find_dmp("nodes.dmp", found)
    if nodes is None:
        return ""
    parent = nodes.parent
    return str(parent)


def attach_taxdump_catalog(payload: Dict[str, Any]) -> Dict[str, Any]:
    """Record each on-disk NCBI dump as ``databases.taxdump.<name:version>``.

    A dump that already has ``nodes.dmp`` does not keep a download recipe.
    ``genomes.taxdump`` stays the default path; ``genomes.taxdump_version`` is
    the stamp of that default.
    """
    genomes = payload.get("genomes") if isinstance(payload.get("genomes"), dict) else {}
    databases = payload.get("databases") if isinstance(payload.get("databases"), dict) else {}
    block = databases.get("taxdump") if isinstance(databases.get("taxdump"), dict) else {}
    rewritten: Dict[str, Any] = {}
    for key, rec in block.items():
        if not isinstance(rec, dict):
            rewritten[str(key)] = rec
            continue
        path = _record_path(rec)
        if path is None or find_dmp("nodes.dmp", path) is None:
            rewritten[str(key)] = rec
            continue
        bare = str(key).split(":", 1)[0] or "ncbi"
        ver = str(rec.get("version") or "").strip() or taxdump_version_label(path)
        item = dict(rec)
        item["path"] = str(path)
        item["type"] = item.get("type") or "database"
        item["lazy-download"] = ""
        if ver:
            item["version"] = ver
        new_key = f"{bare}:{ver}" if ver else str(key)
        rewritten[new_key] = item
    raw = str(genomes.get("taxdump") or "").strip()
    if raw:
        path = _as_path(raw)
        if find_dmp("nodes.dmp", path) is not None:
            ver = str(genomes.get("taxdump_version") or "").strip() or taxdump_version_label(path)
            if ver:
                genomes["taxdump_version"] = ver
            key = f"ncbi:{ver}" if ver else "ncbi"
            rewritten.setdefault(
                key,
                {
                    "path": str(path),
                    "flags": "",
                    "lazy-download": "",
                    "type": "database",
                    "version": ver,
                },
            )
    if rewritten:
        databases = dict(databases)
        databases["taxdump"] = rewritten
        payload["databases"] = databases
    if genomes:
        payload["genomes"] = genomes
    return payload


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
    try:
        from samovar.db_spec import iter_database_records, lookup_database_record

        rec = lookup_database_record(cfg, "taxdump", "ncbi")
        if rec and rec.get("path"):
            return _as_path(rec["path"])
        grouped = iter_database_records(cfg or {}).get("taxdump") or {}
        for key, item in grouped.items():
            bare = str(key).split(":", 1)[0].lower()
            if bare in {"gtdb", "gtdb_taxdump"}:
                continue
            if item.get("path"):
                return _as_path(item["path"])
    except Exception:
        pass
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


def parse_names_dmp(path: PathLike) -> Dict[int, str]:
    """``taxid → scientific name`` from ``names.dmp`` (scientific name rows only)."""
    names: Dict[int, str] = {}
    with open(path, encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            parts = line.split("\t|\t")
            if len(parts) < 4:
                continue
            klass = parts[3].strip().rstrip("|").strip()
            if klass != "scientific name":
                continue
            try:
                taxid = int(parts[0].strip())
            except ValueError:
                continue
            names[taxid] = parts[1].strip()
    return names


def load_scientific_names(path: PathLike) -> Dict[int, str]:
    """Load scientific names with an atomic pickle cache (same idea as nodes.dmp)."""
    names_path = _as_path(path)
    if not names_path.is_file():
        raise FileNotFoundError(f"names.dmp not found: {names_path}")
    try:
        from samovar.parse_annotators import _cache_dir

        cache_dir = _cache_dir()
    except Exception:
        cache_dir = Path(tempfile.gettempdir()) / "samovar_taxdump"
        cache_dir.mkdir(parents=True, exist_ok=True)
    stat = names_path.stat()
    digest = hashlib.sha256(
        f"{names_path.resolve()}:{int(stat.st_mtime)}:{stat.st_size}".encode()
    ).hexdigest()[:16]
    cache_file = cache_dir / f"names_{digest}.pkl"

    def _try_load():
        try:
            with open(cache_file, "rb") as handle:
                return pickle.load(handle)
        except FileNotFoundError:
            return None
        except (OSError, pickle.UnpicklingError, EOFError, AttributeError):
            try:
                cache_file.unlink()
            except OSError:
                pass
            return None

    loaded = _try_load()
    if loaded is not None:
        return loaded
    names = parse_names_dmp(names_path)
    loaded = _try_load()
    if loaded is not None:
        return loaded
    cache_dir.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=cache_file.name + ".", suffix=".tmp", dir=str(cache_dir)
    )
    try:
        with os.fdopen(fd, "wb") as handle:
            pickle.dump(names, handle, protocol=pickle.HIGHEST_PROTOCOL)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_name, cache_file)
    except OSError:
        try:
            os.unlink(tmp_name)
        except OSError:
            pass
        return names
    return names


def load_scientific_name_to_taxid(path: Optional[PathLike] = None) -> Dict[str, str]:
    """Lowercase scientific name → NCBI taxid (smaller id wins on collisions)."""
    names_path = _as_path(path) if path is not None else names_dmp()
    if names_path is None or not Path(names_path).is_file():
        return {}
    names = load_scientific_names(names_path)
    out: Dict[str, str] = {}
    for taxid, name in names.items():
        key = str(name or "").strip().lower()
        if not key:
            continue
        prev = out.get(key)
        if prev is None or int(prev) > int(taxid):
            out[key] = str(taxid)
    return out


_MERGED_CACHE: Optional[Dict[str, str]] = None


def merged_taxid_map(cfg: Optional[Dict[str, Any]] = None) -> Dict[str, str]:
    """``old_taxid → new_taxid`` from ``merged.dmp`` (empty if the file is missing)."""
    global _MERGED_CACHE
    if _MERGED_CACHE is not None:
        return _MERGED_CACHE
    path = find_dmp("merged.dmp", taxdump_dir(cfg))
    mapping: Dict[str, str] = {}
    if path is None:
        _MERGED_CACHE = mapping
        return mapping
    try:
        with path.open(encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                parts = line.split("\t|\t")
                if len(parts) < 2:
                    continue
                old = parts[0].strip()
                new = parts[1].strip().rstrip("|").strip()
                if old.isdigit() and new.isdigit():
                    mapping[old] = new
    except OSError:
        mapping = {}
    _MERGED_CACHE = mapping
    return mapping


def canonical_ncbi_taxid(taxid: str, cfg: Optional[Dict[str, Any]] = None) -> str:
    """Follow ``merged.dmp`` so retired taxids hit the current catalog key."""
    text = str(taxid or "").strip().split(".")[0]
    if not text.isdigit():
        return text
    mapping = merged_taxid_map(cfg)
    seen = set()
    current = text
    while current in mapping and current not in seen:
        seen.add(current)
        current = mapping[current]
    return current


_MERGED: Optional[Dict[str, str]] = None


def merged_map(cfg: Optional[Dict[str, Any]] = None) -> Dict[str, str]:
    """``old_taxid → current_taxid`` from ``merged.dmp`` (cached)."""
    global _MERGED
    if _MERGED is not None:
        return _MERGED
    mapping: Dict[str, str] = {}
    path = find_dmp("merged.dmp", taxdump_dir(cfg))
    if path is None:
        _MERGED = mapping
        return mapping
    try:
        with path.open(encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                parts = [p.strip() for p in line.replace("\t|\n", "").split("\t|")]
                if len(parts) < 2:
                    continue
                old, new = parts[0].strip(), parts[1].strip()
                if old.isdigit() and new.isdigit():
                    mapping[old] = new
    except OSError:
        pass
    _MERGED = mapping
    return mapping


def follow_merged(taxid: str, cfg: Optional[Dict[str, Any]] = None) -> str:
    """Walk ``merged.dmp`` so obsolete taxids hit the current catalog key."""
    text = str(taxid or "").strip().split(".")[0]
    if not text.isdigit():
        return text
    mapping = merged_map(cfg)
    seen = set()
    current = text
    while current in mapping and current not in seen:
        seen.add(current)
        current = mapping[current]
    return current


def taxdump_tarball(cfg: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    roots = [taxdump_dir(cfg)]
    usable = usable_taxdump_dir(cfg)
    if usable is not None:
        roots.insert(0, usable)
    seen = set()
    for root in roots:
        key = str(root)
        if key in seen:
            continue
        seen.add(key)
        for name in ("taxdump.tar.gz", "ncbi-taxonomy.tar.gz"):
            path = root / name
            try:
                if path.is_file():
                    return path
            except OSError:
                continue
        nested = root / "taxonomy"
        for name in ("taxdump.tar.gz", "ncbi-taxonomy.tar.gz"):
            path = nested / name
            try:
                if path.is_file():
                    return path
            except OSError:
                continue
    return None


def local_taxdump_tarball(cfg: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    """Tarball for CAMISIM: the configured archive, or one packed from ``.dmp`` files.

    Packing stays on local disk. This does not fetch NCBI.
    """
    archive = taxdump_tarball(cfg)
    if archive is not None:
        return archive
    root = usable_taxdump_dir(cfg)
    if root is None:
        return None
    dest = root / "taxdump.tar.gz"
    members = []
    for name in TAXDUMP_DMP_FILES:
        found = find_dmp(name, root)
        if found is not None:
            members.append(found)
    if not members:
        return None
    try:
        with tarfile.open(dest, "w:gz") as tar:
            for path in members:
                tar.add(path, arcname=path.name)
    except OSError:
        logger.info("Could not pack a local taxdump tarball in %s", root)
        return None
    return dest if dest.is_file() else None


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
    """Return a directory that contains ``nodes.dmp``.

    A dump already recorded in the config (any version) is reused. NCBI is
    fetched only when none of those directories has ``nodes.dmp``.
    """
    found = usable_taxdump_dir(cfg)
    if found is not None:
        return found
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
