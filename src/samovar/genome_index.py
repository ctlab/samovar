"""Main-config genome index: ``{accession}.fa.gz`` in ``samovar_database``.

``genomes.data`` maps **taxID** to
``[species_level_taxID, genome_ID, database, file_name]``.
Processed files live under ``{samovar_database}/processed/``.
``samovar reindex`` moves ``processed/`` trees into that store and rewrites
the install config pointed at by ``build/config_path``.
"""

from __future__ import annotations

import argparse
import gzip
import logging
import os
import re
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

from samovar.seqio import is_fasta_name, is_gzip_path, sequence_stem

logger = logging.getLogger(__name__)

PathLike = Union[str, os.PathLike]

ACCESSION_RE = re.compile(r"^(GC[AF]_\d+(?:\.\d+)?)$", re.IGNORECASE)
FOLDER_ID = "samovar_database"
TEST_FOLDER_ID = "test"
TEST_TAXID_SUFFIX = "_test"
RECORD_LEN = 4  # species_taxid, genome_id, database, file_name


def _as_path(path: PathLike) -> Path:
    return Path(path).expanduser()


def is_test_taxid(value: str) -> bool:
    """True for catalog keys like ``562_test`` (bundled ISS stubs, not NCBI)."""
    text = str(value or "").strip()
    if not text.endswith(TEST_TAXID_SUFFIX):
        return False
    return text[: -len(TEST_TAXID_SUFFIX)].split(".")[0].isdigit()


def numeric_taxid_stem(value: str) -> str:
    text = str(value or "").strip()
    if is_test_taxid(text):
        text = text[: -len(TEST_TAXID_SUFFIX)]
    return text.split(".")[0]


def as_test_taxid(value: str) -> str:
    stem = numeric_taxid_stem(value)
    if stem.isdigit():
        return f"{stem}{TEST_TAXID_SUFFIX}"
    text = str(value or "").strip()
    if text and not text.endswith(TEST_TAXID_SUFFIX):
        return f"{text}{TEST_TAXID_SUFFIX}"
    return text


def _is_taxid(value: str) -> bool:
    text = str(value or "").strip()
    if not text or is_assembly_accession(text):
        return False
    if is_test_taxid(text):
        return True
    return text.split(".")[0].isdigit()


def _is_ncbi_taxid(value: str) -> bool:
    """Numeric NCBI taxid (not ``562_test``)."""
    text = str(value or "").strip()
    if not text or is_test_taxid(text) or is_assembly_accession(text):
        return False
    return text.split(".")[0].isdigit()


def _looks_like_fasta_name(name: str) -> bool:
    lower = str(name or "").lower()
    return any(
        lower.endswith(ext)
        for ext in (
            ".fa",
            ".fa.gz",
            ".fna",
            ".fna.gz",
            ".fasta",
            ".fasta.gz",
            ".fsa",
            ".fsa.gz",
        )
    )


def coerce_genome_record(key: str, rec: Sequence[Any]) -> Tuple[str, List[str]]:
    """Return ``(taxid_key, [species, genome_id, database, file_name])``.

    Accepts the current 4-field schema and older catalogs:
    ``[accession, folder_id, file_name]`` and
    ``[accession, folder_id, file_name, taxid]``.
    """
    key = str(key or "").strip()
    parts = [str(x).strip() if x is not None else "" for x in rec]
    species = genome_id = database = file_name = taxid = ""

    if len(parts) >= 4 and not _looks_like_fasta_name(parts[2]):
        # Current: [species_taxid, genome_id, database, file_name]
        species, genome_id, database, file_name = parts[0], parts[1], parts[2], parts[3]
        taxid = key if _is_taxid(key) else (species if _is_taxid(species) else "")
    elif len(parts) >= 4 and _is_taxid(parts[3]) and (
        is_assembly_accession(parts[0]) or is_assembly_accession(key)
    ):
        # Intermediate: [accession, folder_id, file_name, taxid]
        genome_id, database, file_name, taxid = parts[0], parts[1], parts[2], parts[3]
        species = taxid
    elif len(parts) >= 3:
        # Legacy: [accession_or_empty, folder_id, file_name]
        genome_id, database, file_name = parts[0], parts[1], parts[2]
        if _is_taxid(key):
            taxid = key
            species = key
            if not genome_id:
                genome_id = key
        elif is_assembly_accession(key) or is_assembly_accession(genome_id):
            genome_id = genome_id or key
        else:
            taxid = key
            species = key
            if not genome_id:
                genome_id = key
    elif parts:
        genome_id = parts[0]
        database = FOLDER_ID
        file_name = processed_filename(key or genome_id)

    if not taxid:
        if _is_taxid(key):
            taxid = key
        elif _is_taxid(species):
            taxid = species
        else:
            taxid = genome_id or key
    if not species:
        species = taxid if _is_taxid(taxid) else ""
    if not file_name:
        file_name = processed_filename(genome_id or taxid or key)
    return taxid, [species, genome_id, database, file_name]


def _record_is_test_stub(key: str, row: Sequence[str]) -> bool:
    species, genome_id, database, file_name = (list(row) + [""] * 4)[:4]
    if database == TEST_FOLDER_ID or is_test_taxid(key):
        return True
    blob = f"{file_name} {genome_id} {species}"
    return "test_genomes" in blob.replace("\\", "/")


def normalize_genome_data(raw: Any) -> Dict[str, List[str]]:
    """Rewrite ``genomes.data`` to taxID → 4-field records.

    Bundled ``data/test_genomes`` stubs are stored as ``{taxid}_test`` so they
    never collide with a real NCBI taxid (and are not used as a library hit).
    """
    if not isinstance(raw, dict):
        return {}
    out: Dict[str, List[str]] = {}
    for key, rec in raw.items():
        if str(key).startswith("_"):
            continue
        if isinstance(rec, (list, tuple)):
            taxid, row = coerce_genome_record(str(key), rec)
        elif rec:
            taxid, row = coerce_genome_record(str(key), [rec, FOLDER_ID, processed_filename(str(key))])
        else:
            continue
        if not taxid:
            continue
        if _record_is_test_stub(taxid, row) or _record_is_test_stub(str(key), row):
            marked = as_test_taxid(taxid if _is_ncbi_taxid(taxid) or is_test_taxid(taxid) else str(key))
            species, genome_id, database, file_name = (row + [""] * 4)[:4]
            if _is_ncbi_taxid(species) or not species:
                species = marked
            if _is_ncbi_taxid(genome_id) or not genome_id:
                genome_id = marked
            row = [species, genome_id, TEST_FOLDER_ID, file_name]
            taxid = marked
        out[taxid] = row
    return out


def normalize_accession(value: str) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    name = Path(text).name
    stem = sequence_stem(name) if "." in name and not ACCESSION_RE.match(name) else name
    if ACCESSION_RE.match(stem):
        return stem
    if ACCESSION_RE.match(name):
        return name
    return text


def is_assembly_accession(value: str) -> bool:
    return bool(ACCESSION_RE.match(str(value or "").strip()))


def processed_filename(accession: str) -> str:
    acc = normalize_accession(accession) or str(accession).strip()
    return f"{acc}.fa.gz"


def raw_filename(accession: str) -> str:
    acc = normalize_accession(accession) or str(accession).strip()
    return f"{acc}.fna.gz"


def _cfg():
    from samovar.paths import load_config

    return load_config()


def samovar_database_dir(cfg: Optional[Dict[str, Any]] = None) -> Path:
    """Install processed-genome store (``genomes.samovar_database``)."""
    from samovar.main_config import genomes_block
    from samovar.paths import repo_root

    data = cfg if cfg is not None else _cfg()
    block = genomes_block(data) if data else {}
    raw = str(
        (block.get("samovar_database") if isinstance(block, dict) else "")
        or data.get("samovar_database")
        or ""
    ).strip()
    if raw:
        return Path(raw).expanduser()
    env = os.environ.get("SAMOVAR_DATABASE", "").strip()
    if env:
        return Path(env).expanduser()
    return repo_root() / "genomes"


def processed_storage_dir(cfg: Optional[Dict[str, Any]] = None) -> Path:
    return samovar_database_dir(cfg) / "processed"


def proteome_storage_dir(cfg: Optional[Dict[str, Any]] = None) -> Path:
    return samovar_database_dir(cfg) / "proteomes"


def raw_storage_dir(cfg: Optional[Dict[str, Any]] = None) -> Path:
    return samovar_database_dir(cfg) / "raw"


def run_genomes_dir(output_dir: PathLike) -> Path:
    return _as_path(output_dir) / ".genomes"


def run_processed_dir(output_dir: PathLike) -> Path:
    return run_genomes_dir(output_dir) / "processed"


def run_raw_dir(output_dir: PathLike) -> Path:
    return run_genomes_dir(output_dir) / "raw"


def folder_id_path(folder_id: str, cfg: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    from samovar.main_config import folder_map, genomes_block

    data = cfg if cfg is not None else _cfg()
    block = genomes_block(data)
    fid = str(folder_id or "").strip()
    if fid in {"samovar_database", "default", "processed"}:
        return processed_storage_dir(data)
    if fid == "raw":
        return raw_storage_dir(data)
    mapped = folder_map(block.get("processed"))
    mapped.update(folder_map(block.get("raw")))
    path = mapped.get(fid)
    if path:
        return Path(path).expanduser()
    tests = block.get("test") or []
    if fid == "test" and tests:
        return Path(str(tests[0])).expanduser()
    return None


def genome_data_map(cfg: Optional[Dict[str, Any]] = None) -> Dict[str, List[str]]:
    from samovar.main_config import genomes_block

    data = cfg if cfg is not None else _cfg()
    block = genomes_block(data)
    return normalize_genome_data(block.get("data") if isinstance(block.get("data"), dict) else {})


def indexed_record(ident: str, cfg: Optional[Dict[str, Any]] = None) -> Optional[Tuple[str, List[str]]]:
    ident = str(ident or "").strip()
    if not ident:
        return None
    data = genome_data_map(cfg)
    acc = normalize_accession(ident)
    if ident in data:
        return ident, data[ident]
    if acc and acc in data:
        return acc, data[acc]
    stem = sequence_stem(ident)
    if stem in data:
        return stem, data[stem]
    want_test = is_test_taxid(ident)
    for key, rec in data.items():
        if is_test_taxid(key) and not want_test:
            # Stubs must not satisfy a real NCBI taxid / accession lookup.
            continue
        species, genome_id, _database, file_name = (rec + [""] * 4)[:4]
        if genome_id and genome_id in {ident, acc, stem}:
            return key, rec
        if species and species == ident:
            return key, rec
        if file_name and sequence_stem(file_name) in {ident, acc, stem}:
            return key, rec
    return None


def resolve_indexed_path(ident: str, cfg: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    hit = indexed_record(ident, cfg)
    if hit is None:
        return None
    _key, rec = hit
    _species, genome_id, database, file_name = (rec + [""] * 4)[:4]
    name = file_name or processed_filename(genome_id or ident)
    base = folder_id_path(database, cfg)
    candidates: List[Path] = []
    if base is not None:
        candidates.append(base / name)
        if base.name != "processed":
            candidates.append(base / "processed" / name)
    candidates.append(processed_storage_dir(cfg) / name)
    for path in candidates:
        try:
            if path.is_file():
                return path
        except OSError:
            continue
    return None


def ensure_gzip_fasta(path: Path) -> Path:
    """Return a ``.gz`` FASTA; gzip in place if the file is still plain.

    Bundled ``data/test_genomes`` files are never rewritten (no gzip, no move).
    """
    src = _as_path(path)
    if not src.is_file():
        raise FileNotFoundError(src)
    from samovar.genome_cache import is_bundled_test_genomes_path

    if is_bundled_test_genomes_path(src):
        return src
    if is_gzip_path(src):
        return src
    dest = src.with_name(src.name + ".gz")
    with src.open("rb") as fin, gzip.open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    src.unlink()
    return dest


def iter_processed_fastas(root: PathLike) -> List[Path]:
    """FASTAs under ``processed/`` (or ``.genomes/processed``) of ``root``."""
    base = _as_path(root)
    if not base.exists():
        return []
    search: List[Path] = []
    if base.is_file():
        return [base] if is_fasta_name(base.name, protein=False) else []
    named = [
        base / "processed",
        base / ".genomes" / "processed",
        base / "genomes" / "processed",
    ]
    if base.name == "processed":
        named.insert(0, base)
    found: List[Path] = []
    seen = set()

    def add(path: Path) -> None:
        try:
            key = str(path.resolve())
        except OSError:
            key = str(path)
        if key in seen:
            return
        if not path.is_file():
            return
        if not is_fasta_name(path.name, protein=False, nucleotide=True):
            return
        seen.add(key)
        found.append(path)

    for directory in named:
        if directory.is_dir():
            search.append(directory)
    if not search and base.is_dir():
        # Directory of processed FASTAs with no subdirectory.
        search.append(base)
    for directory in search:
        try:
            children = list(directory.iterdir())
        except OSError:
            continue
        for child in sorted(children):
            add(child)
    return found


def accession_from_fasta(path: PathLike) -> str:
    name = _as_path(path).name
    stem = sequence_stem(name)
    if stem.endswith("-processed"):
        stem = stem[: -len("-processed")]
    acc = normalize_accession(stem)
    return acc or stem


def _write_data_map(data: Dict[str, List[str]]) -> None:
    from samovar.main_config import genomes_block
    from samovar.paths import load_config, update_config

    cfg = load_config()
    block = genomes_block(cfg) or {}
    block["data"] = data
    store = str(samovar_database_dir(cfg))
    block["samovar_database"] = store
    processed = dict(block.get("processed") or {})
    processed[FOLDER_ID] = str(processed_storage_dir(cfg))
    processed["default"] = processed.get("default") or str(processed_storage_dir(cfg))
    block["processed"] = processed
    update_config({"genomes": block})


def stage_into_dir(src: PathLike, dest_dir: PathLike) -> Path:
    """Hardlink (fallback: copy) a processed FASTA into ``dest_dir`` under ``{acc}.fa.gz``."""
    src_path = ensure_gzip_fasta(_as_path(src))
    dest = _as_path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    target = dest / processed_filename(accession_from_fasta(src_path))
    try:
        same = src_path.resolve() == target.resolve()
    except OSError:
        same = str(src_path) == str(target)
    if same:
        return target
    if target.exists() or target.is_symlink():
        target.unlink()
    try:
        os.link(src_path, target)
    except OSError:
        shutil.copy2(src_path, target)
    return target


def numeric_taxid_for(ident: str, cfg: Optional[Dict[str, Any]] = None) -> str:
    """NCBI taxonomy id for a taxid or assembly accession recorded in the index."""
    ident = str(ident or "").strip()
    if not ident:
        return ""
    if is_test_taxid(ident):
        return ident
    if _is_ncbi_taxid(ident):
        return ident.split(".")[0]
    hit = indexed_record(ident, cfg)
    if hit is not None:
        key, rec = hit
        if is_test_taxid(key):
            return key
        if _is_ncbi_taxid(key):
            return key.split(".")[0]
        species = rec[0] if rec else ""
        if is_test_taxid(species):
            return species
        if _is_ncbi_taxid(species):
            return species.split(".")[0]
    return ident


def index_processed_file(
    path: PathLike,
    *,
    accession: str = "",
    folder_id: str = "",
    move_to: Optional[Path] = None,
    taxid: str = "",
    species_taxid: str = "",
) -> Path:
    """Gzip if needed, optionally move into ``move_to``, and record in config."""
    from samovar.main_config import folder_map, genomes_block
    from samovar.paths import load_config, update_config

    src = _as_path(path)
    from samovar.genome_cache import is_bundled_test_genomes_path

    if is_bundled_test_genomes_path(src):
        try:
            resolved = src.resolve()
        except OSError:
            resolved = src
        rel = src.name
        parts = resolved.parts
        if "test_genomes" in parts:
            idx = parts.index("test_genomes")
            rel = "/".join(parts[idx + 1 :]) or src.name
        else:
            try:
                from samovar.paths import test_genomes_dir as _tgd

                rel = resolved.relative_to(_tgd().resolve()).as_posix()
            except (OSError, ValueError):
                rel = src.name
        taxid_key = as_test_taxid(
            str(taxid or species_taxid or accession or accession_from_fasta(src) or src.stem)
        )
        cfg = load_config()
        block = genomes_block(cfg) or {}
        data = genome_data_map(cfg)
        for stale in [k for k in list(data) if k in {taxid_key, numeric_taxid_stem(taxid_key)}]:
            if _record_is_test_stub(stale, data.get(stale) or []):
                data.pop(stale, None)
        data[taxid_key] = [taxid_key, taxid_key, TEST_FOLDER_ID, rel]
        block["data"] = data
        update_config({"genomes": block})
        logger.info("Recorded test genome %s -> %s (not copied into the library)", taxid_key, src)
        return src

    src = ensure_gzip_fasta(src)
    acc = accession or accession_from_fasta(src)
    dest_dir = _as_path(move_to) if move_to is not None else src.parent
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / processed_filename(acc)
    if src.resolve() != dest.resolve():
        if dest.exists():
            dest.unlink()
        shutil.move(str(src), str(dest))
        src = dest
    elif src.name != dest.name:
        src = src.rename(dest)
    store = processed_storage_dir()
    try:
        in_store = dest_dir.resolve() == store.resolve()
    except OSError:
        in_store = False
    fid = folder_id or (FOLDER_ID if in_store else dest_dir.name or "lib")
    reserved = {FOLDER_ID, "default", "processed", "raw", "test"}
    if not in_store and fid in reserved:
        fid = dest_dir.parent.name or "lib"
        if fid in reserved:
            fid = f"lib_{abs(hash(str(dest_dir.resolve())))}"
    cfg = load_config()
    block = genomes_block(cfg) or {}
    proc = folder_map(block.get("processed"))
    proc[fid] = str(dest_dir)
    proc.setdefault("default", str(store))
    proc.setdefault(FOLDER_ID, str(store))
    block["processed"] = proc
    block["samovar_database"] = str(samovar_database_dir(cfg))
    data = genome_data_map(cfg)
    prev = indexed_record(acc) or indexed_record(taxid)
    prev_row = prev[1] if prev else []
    organism = str(taxid or "").strip()
    species = str(species_taxid or "").strip()
    if not organism and prev_row:
        organism = prev[0] if prev and _is_taxid(prev[0]) else ""
        if not organism and _is_taxid(prev_row[0]):
            organism = prev_row[0]
    if not species:
        species = (prev_row[0] if prev_row and _is_taxid(prev_row[0]) else "") or organism
    taxid_key = organism if _is_taxid(organism) else (acc or dest.stem)
    stale = [
        key
        for key, row in data.items()
        if key in {acc, taxid_key, organism}
        or (row and row[1] == acc)
    ]
    for key in stale:
        data.pop(key, None)
    data[taxid_key] = [species or organism, acc, fid, dest.name]
    block["data"] = data
    update_config({"genomes": block})
    logger.info("Indexed %s (%s) -> %s", taxid_key, acc, src)
    return src


def catalog_enabled(cfg: Optional[Dict[str, Any]] = None) -> bool:
    """True when the install config (or ``SAMOVAR_DATABASE``) names a genome store."""
    from samovar.main_config import genomes_block

    data = cfg if cfg is not None else _cfg()
    block = genomes_block(data) if data else {}
    if str((block or {}).get("samovar_database") or "").strip():
        return True
    return bool(os.environ.get("SAMOVAR_DATABASE", "").strip())


def catalog_processed_genome(
    path: PathLike,
    *,
    taxid: str = "",
    accession: str = "",
    species_taxid: str = "",
    keep_src: bool = True,
    stage_dir: Optional[PathLike] = None,
) -> Path:
    """Copy/hardlink ``path`` into ``samovar_database/processed`` and update config.

    Run-folder FASTAs are left in place when ``keep_src`` is true; the catalog
    always points at the library copy so later Kaiju/Kraken/ISS stages reuse it.
    """
    if not catalog_enabled():
        src = _as_path(path)
        if stage_dir:
            return stage_into_dir(src, stage_dir)
        return src
    src = ensure_gzip_fasta(_as_path(path))
    acc = accession or accession_from_fasta(src)
    store = processed_storage_dir()
    store.mkdir(parents=True, exist_ok=True)
    target = store / processed_filename(acc)
    try:
        same = src.resolve() == target.resolve()
    except OSError:
        same = str(src) == str(target)
    if not same:
        if target.exists() or target.is_symlink():
            try:
                if target.stat().st_ino == src.stat().st_ino:
                    same = True
            except OSError:
                same = False
        if not same:
            if target.exists() or target.is_symlink():
                target.unlink()
            if keep_src:
                try:
                    os.link(src, target)
                except OSError:
                    shutil.copy2(src, target)
            else:
                shutil.move(str(src), str(target))
            src = target
    indexed = index_processed_file(
        src,
        accession=acc,
        taxid=taxid,
        species_taxid=species_taxid,
        move_to=store,
    )
    if stage_dir:
        return stage_into_dir(indexed, stage_dir)
    return indexed


def harvest_run_genomes(output_dir: PathLike) -> int:
    """Index processed FASTAs under a run dir into the install catalog."""
    if not catalog_enabled():
        return 0
    counted = 0
    roots = [
        _as_path(output_dir) / ".genomes",
    ]
    seen = set()
    for root in roots:
        for path in iter_processed_fastas(root):
            try:
                key = str(path.resolve())
            except OSError:
                key = str(path)
            if key in seen:
                continue
            seen.add(key)
            try:
                catalog_processed_genome(path, keep_src=True, stage_dir=path.parent)
                counted += 1
            except Exception:
                logger.warning("Could not catalog %s", path, exc_info=True)
    if counted:
        logger.info("Harvested %s processed genome(s) into the install catalog", counted)
    return counted


def register_database(tool: str, name: str, path: str, flags: str = "") -> None:
    from samovar.main_config import databases_of
    from samovar.paths import load_config, update_config

    cfg = load_config()
    dbs = databases_of(cfg)
    rows = list(dbs.get(tool) or [])
    entry = [str(name).strip(), str(path).strip(), str(flags or "").strip()]
    replaced = False
    for i, row in enumerate(rows):
        if row and row[0] == entry[0]:
            rows[i] = entry
            replaced = True
            break
    if not replaced:
        rows.append(entry)
    dbs[str(tool)] = rows
    update_config({"databases": dbs})
    logger.info("Registered database %s/%s -> %s", tool, name, path)


def lookup_database(tool: str, name: str) -> Optional[List[str]]:
    from samovar.main_config import databases_of
    from samovar.paths import load_config

    rows = databases_of(load_config()).get(tool) or []
    for row in rows:
        if row and row[0] == name:
            return row
    return None


def reindex(
    sources: Optional[Sequence[PathLike]] = None,
    dest: Optional[PathLike] = None,
) -> Dict[str, Any]:
    """Move processed genomes into ``dest`` (default samovar_database) and update config.

    * With ``sources``: harvest ``processed/`` trees (error if none).
    * With no sources: move every already-indexed FASTA that is not yet in dest.
      Error if the index has no processed files on disk.
    """
    dest_dir = _as_path(dest) if dest else processed_storage_dir()
    dest_dir.mkdir(parents=True, exist_ok=True)
    moved: List[str] = []
    from samovar.genome_cache import is_bundled_test_genomes_path

    if sources:
        files: List[Path] = []
        for src in sources:
            files.extend(iter_processed_fastas(src))
        files = [path for path in files if not is_bundled_test_genomes_path(path)]
        if not files:
            raise RuntimeError(
                "reindex: nothing in processed under "
                + ", ".join(str(s) for s in sources)
            )
        for path in files:
            acc = accession_from_fasta(path)
            dest_file = dest_dir / processed_filename(acc)
            try:
                if dest_file.is_file() and dest_file.stat().st_ino == path.stat().st_ino:
                    continue
            except OSError:
                pass
            indexed = resolve_indexed_path(acc)
            if indexed is not None:
                try:
                    if indexed.resolve().parent == dest_dir.resolve() and is_gzip_path(indexed):
                        continue
                except OSError:
                    pass
            placed = index_processed_file(path, move_to=dest_dir)
            moved.append(str(placed))
        return {"dest": str(dest_dir), "moved": moved, "count": len(moved)}

    data = genome_data_map()
    files = []
    for key, rec in data.items():
        if is_test_taxid(key) or _record_is_test_stub(key, rec):
            continue
        path = resolve_indexed_path(key)
        if path is None or is_bundled_test_genomes_path(path):
            continue
        try:
            if path.resolve().parent == dest_dir.resolve() and is_gzip_path(path):
                continue
        except OSError:
            pass
        files.append(path)
    if not files:
        extra = iter_processed_fastas(processed_storage_dir())
        if not extra:
            raise RuntimeError("reindex: nothing in processed")
        return {"dest": str(dest_dir), "moved": [], "count": 0, "already": len(extra)}
    for path in files:
        placed = index_processed_file(path, move_to=dest_dir)
        moved.append(str(placed))
    return {"dest": str(dest_dir), "moved": moved, "count": len(moved)}


def drop_indexed(*idents: str) -> None:
    data = genome_data_map()
    changed = False
    for ident in idents:
        rec = indexed_record(ident)
        if rec is None:
            continue
        key, _ = rec
        data.pop(key, None)
        changed = True
    if changed:
        _write_data_map(data)


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Move processed genomes into samovar_database and update the install config"
    )
    parser.add_argument(
        "sources",
        nargs="*",
        help="Folders with genomes/processed (default: already-indexed genomes)",
    )
    parser.add_argument(
        "--dest",
        default=None,
        help="Processed storage (default: genomes.samovar_database/processed)",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    args = parse_args(argv)
    try:
        result = reindex(args.sources or None, dest=args.dest)
    except RuntimeError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(f"reindex dest={result['dest']} moved={result['count']}")
    for path in result.get("moved") or []:
        print(f"  {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
