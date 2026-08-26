"""Main-config genome index: ``{accession}.fa.gz`` in ``samovar_database``.

``genomes.data`` maps accession (preferred) or taxid to
``[accession, folder_id, file_name]``. Processed files live under
``{samovar_database}/processed/``. ``samovar reindex`` moves ``processed/``
trees into that store and rewrites the install config pointed at by
``build/config_path``.
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


def _as_path(path: PathLike) -> Path:
    return Path(path).expanduser()


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
    raw = block.get("data") if isinstance(block.get("data"), dict) else {}
    out: Dict[str, List[str]] = {}
    for key, rec in raw.items():
        if str(key).startswith("_"):
            continue
        if isinstance(rec, (list, tuple)):
            parts = [str(x).strip() if x is not None else "" for x in rec]
            while len(parts) < 3:
                parts.append("")
            out[str(key)] = parts[:4] if len(parts) > 3 else parts[:3]
        elif rec:
            out[str(key)] = [str(rec).strip(), FOLDER_ID, processed_filename(str(key))]
    return out


def indexed_record(ident: str, cfg: Optional[Dict[str, Any]] = None) -> Optional[Tuple[str, List[str]]]:
    ident = str(ident or "").strip()
    if not ident:
        return None
    data = genome_data_map(cfg)
    acc = normalize_accession(ident)
    if acc and acc in data:
        return acc, data[acc]
    if ident in data:
        return ident, data[ident]
    stem = sequence_stem(ident)
    if stem in data:
        return stem, data[stem]
    for key, rec in data.items():
        if rec and rec[0] and rec[0] == acc:
            return key, rec
        if rec and rec[2] and sequence_stem(rec[2]) in {ident, acc, stem}:
            return key, rec
    return None


def resolve_indexed_path(ident: str, cfg: Optional[Dict[str, Any]] = None) -> Optional[Path]:
    hit = indexed_record(ident, cfg)
    if hit is None:
        return None
    _key, rec = hit
    accession, folder_id, file_name = rec[0], rec[1], rec[2]
    name = file_name or processed_filename(accession or ident)
    base = folder_id_path(folder_id, cfg)
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
    """Return a ``.gz`` FASTA; gzip in place if the file is still plain."""
    src = _as_path(path)
    if not src.is_file():
        raise FileNotFoundError(src)
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
    if ident.split(".")[0].isdigit() and not is_assembly_accession(ident):
        return ident.split(".")[0]
    hit = indexed_record(ident, cfg)
    if hit is not None:
        rec = hit[1]
        if len(rec) > 3 and str(rec[3]).split(".")[0].isdigit():
            return str(rec[3]).split(".")[0]
    return ident


def index_processed_file(
    path: PathLike,
    *,
    accession: str = "",
    folder_id: str = "",
    move_to: Optional[Path] = None,
    taxid: str = "",
) -> Path:
    """Gzip if needed, optionally move into ``move_to``, and record in config."""
    from samovar.main_config import folder_map, genomes_block
    from samovar.paths import load_config, update_config

    src = ensure_gzip_fasta(_as_path(path))
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
    row = [acc, fid, dest.name]
    if taxid:
        row.append(str(taxid).strip())
    elif acc in data and len(data[acc]) > 3:
        row.append(data[acc][3])
    data[acc] = row
    block["data"] = data
    update_config({"genomes": block})
    logger.info("Indexed %s -> %s", acc, src)
    return src


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
    if sources:
        files: List[Path] = []
        for src in sources:
            files.extend(iter_processed_fastas(src))
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
        path = resolve_indexed_path(key)
        if path is None:
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
