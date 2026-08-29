"""
Genome fetching and taxonomy parsing functionality
"""

import os
import logging
import shutil
import socket
from typing import Any, Collection, FrozenSet, List, Optional, Sequence, Tuple
import urllib.request
from Bio import Entrez
from pathlib import Path
import random
from tqdm import tqdm
import time
from samovar.fasta_processor import preprocess_fasta
from samovar.genome_cache import (
    allow_test_genomes,
    find_library_genome,
    genome_download_dir,
    place_genome,
    processed_genomes_dir as cache_processed_dir,
    register_genome_dir,
    reuse_enabled,
    ensure_genome_config,
)
from samovar.seqio import (
    find_existing_processed_genome,
    gunzip_file,
    gzip_file,
    is_gzip_path,
    processed_genome_path,
)
from samovar.genome_index import (
    accession_from_fasta,
    catalog_processed_genome,
    is_assembly_accession,
    normalize_accession,
    processed_filename,
    processed_storage_dir,
    raw_filename,
    raw_storage_dir,
    resolve_indexed_path,
    run_processed_dir,
    run_raw_dir,
    index_processed_file,
    proteome_storage_dir,
    stage_into_dir,
)

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def default_entrez_email() -> str:
    """NCBI Entrez contact email from the environment or install config."""
    from samovar.paths import ncbi_email

    return ncbi_email()

def _entrez_retry(func, max_retries=3, initial_delay=1):
    """
    Retry an Entrez function with exponential backoff.
    
    Args:
        func: Function to retry
        max_retries: Maximum number of retries
        initial_delay: Initial delay in seconds
        
    Returns:
        Result of the function call
    """
    delay = initial_delay
    last_exception = None
    
    for attempt in range(max_retries):
        try:
            return func()
        except Exception as e:
            last_exception = e
            msg = str(e)
            rate_limited = (
                "429" in msg
                or "Too Many Requests" in msg
                or "HTTP Error 400" in msg  # NCBI sometimes returns 400 under load
                or "Remote end closed connection" in msg
            )
            if rate_limited and attempt < max_retries - 1:
                logger.warning(f"Rate limited / transient NCBI error, retrying in {delay} seconds...")
                time.sleep(delay)
                delay = min(delay * 2, 8)
                continue
            raise
    
    raise last_exception

def bundled_test_genome(taxid: str | int) -> Optional[Path]:
    """Return a genome shipped under ``data/test_genomes`` for this taxid, if any.

    These files are truncated ISS/CI stubs. ``fetch_genome`` does not use them
    unless ``--test-genomes`` / ``SAMOVAR_ALLOW_TEST_GENOMES=1`` is set.
    """
    from samovar.paths import test_genomes_dir

    taxid = _normalize_taxid(taxid)
    root = test_genomes_dir()
    if not root.is_dir():
        return None
    names = [
        f"{taxid}{ext}"
        for ext in (".fna", ".fa", ".fasta", ".fna.gz", ".fa.gz", ".fasta.gz")
    ]
    search_roots = [root, root / "meta", root / "host"]
    for base in search_roots:
        if not base.is_dir():
            continue
        for name in names:
            candidate = base / name
            if candidate.is_file():
                return candidate
    return None


def _copy_bundled_genome(taxid: str | int, output_folder: str) -> Optional[str]:
    bundled = bundled_test_genome(taxid)
    if bundled is None:
        return None
    dest_dir = Path(output_folder)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest = dest_dir / f"{_normalize_taxid(taxid)}{''.join(bundled.suffixes)}"
    if not dest.exists():
        shutil.copy2(bundled, dest)
    return str(dest)


def _download_url(url: str, dest: Path, timeout: int = 45) -> None:
    """Fetch ``url`` to ``dest`` with a socket timeout (urlretrieve has none)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    req = urllib.request.Request(url)
    expected: Optional[int] = None
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as out:
        cl = resp.headers.get("Content-Length")
        if cl:
            try:
                expected = int(cl)
            except (TypeError, ValueError):
                expected = None
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    if expected is not None:
        warn_gzip_size_vs_ncbi(dest, identity=dest.name, expected_bytes=expected)


def _normalize_taxid(taxid: str | int) -> str:
    taxid = str(taxid)
    return taxid.split(".")[0]


UNLIMITED_GENOME_MB = float("inf")


def parse_max_genome_mb(
    value: Any = None,
    *,
    default_from_env: bool = True,
) -> float:
    """Size cap for a *new* NCBI download. ``0`` / ``inf`` / empty → unlimited."""
    if value is None and default_from_env:
        value = os.environ.get("SAMOVAR_MAX_GENOME_MB", "")
    if value is None or value is False:
        return UNLIMITED_GENOME_MB
    if isinstance(value, str):
        text = value.strip().lower()
        if text in {"", "inf", "infinity", "+inf", "none", "unlimited", "unlimited.0"}:
            return UNLIMITED_GENOME_MB
        try:
            value = float(text)
        except ValueError:
            return UNLIMITED_GENOME_MB
    try:
        number = float(value)
    except (TypeError, ValueError):
        return UNLIMITED_GENOME_MB
    if number <= 0 or number == UNLIMITED_GENOME_MB:
        return UNLIMITED_GENOME_MB
    return number


def argparse_max_genome_mb(text: str) -> float:
    return parse_max_genome_mb(text, default_from_env=False)


def parse_genome_skip_list(
    value: Any = None,
    *,
    default_from_env: bool = True,
) -> FrozenSet[str]:
    """Taxids and/or assembly accessions that must not be newly downloaded."""
    if value is None and default_from_env:
        value = os.environ.get("SAMOVAR_GENOME_SKIP_LIST", "")
    if value is None or value is False:
        return frozenset()
    if isinstance(value, (list, tuple, set, frozenset)):
        parts = [str(item) for item in value]
    else:
        text = str(value).strip()
        if not text:
            return frozenset()
        parts = [
            piece
            for piece in text.replace(";", ",").replace(":", ",").replace("\n", ",").split(",")
            if piece.strip()
        ]
    keys: set[str] = set()
    for part in parts:
        token = str(part).strip()
        if not token:
            continue
        keys.add(token)
        keys.add(_normalize_taxid(token))
        acc = normalize_accession(token)
        if acc:
            keys.add(acc)
    keys.discard("")
    return frozenset(keys)


def genome_id_in_skip_list(identity: str, skip: Collection[str]) -> bool:
    if not skip:
        return False
    raw = str(identity or "").strip()
    if not raw:
        return False
    candidates = {raw, _normalize_taxid(raw)}
    acc = normalize_accession(raw)
    if acc:
        candidates.add(acc)
    return bool(candidates & set(skip))


def _genomic_fna_gz_url(ftp_path: str) -> str:
    asm_name = os.path.basename(str(ftp_path).rstrip("/"))
    return f"{ftp_path}/{asm_name}_genomic.fna.gz".replace("ftp://", "https://")


def _remote_content_length_bytes(url: str, timeout: int = 30) -> Optional[int]:
    if not url:
        return None
    try:
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            cl = resp.headers.get("Content-Length")
            if cl:
                return int(cl)
    except Exception:
        return None
    return None


def _remote_fasta_size_bytes(ftp_path: str, timeout: int = 30) -> Optional[int]:
    """HTTP HEAD Content-Length of ``{asm}_genomic.fna.gz``, or None if unknown."""
    if not ftp_path:
        return None
    return _remote_content_length_bytes(_genomic_fna_gz_url(ftp_path), timeout=timeout)


def _remote_fasta_size_mb(ftp_path: str, timeout: int = 30) -> Optional[float]:
    nbytes = _remote_fasta_size_bytes(ftp_path, timeout=timeout)
    if nbytes is None:
        return None
    return nbytes / (1024 * 1024)


_GZIP_SIZE_REL_TOL = 0.01
_GZIP_SIZE_ABS_TOL = 64


def gzip_sizes_disagree(local_bytes: int, expected_bytes: int) -> bool:
    """True when on-disk gzip size is not close to NCBI Content-Length."""
    if expected_bytes <= 0:
        return local_bytes != expected_bytes
    delta = abs(int(local_bytes) - int(expected_bytes))
    return delta > max(_GZIP_SIZE_ABS_TOL, int(expected_bytes * _GZIP_SIZE_REL_TOL))


def warn_gzip_size_vs_ncbi(
    path: Path,
    *,
    identity: str = "",
    ftp_path: Optional[str] = None,
    url: Optional[str] = None,
    expected_bytes: Optional[int] = None,
    timeout: int = 30,
) -> None:
    """Warn if a saved ``.gz`` is much smaller/larger than NCBI Content-Length.

    Existing cache files are still reused; this is integrity signal only.
    """
    target = Path(path)
    if not target.is_file():
        return
    local = target.stat().st_size
    expected = expected_bytes
    if expected is None and ftp_path:
        expected = _remote_fasta_size_bytes(ftp_path, timeout=timeout)
    if expected is None and url:
        expected = _remote_content_length_bytes(url, timeout=timeout)
    if expected is None:
        return
    if gzip_sizes_disagree(local, expected):
        logger.warning(
            "Cached gzip for %s is %s bytes on disk but NCBI reports %s bytes "
            "(%s); the file may be truncated, incomplete, or stale",
            identity or target.name,
            local,
            expected,
            target,
        )


def _warn_existing_ncbi_gzip(path: Path, taxid: str, email: str) -> None:
    if not is_gzip_path(path):
        return
    name = path.name.lower()
    if not (name.endswith(".fna.gz") or name.endswith(".fasta.gz")):
        return
    try:
        ftp_path = _assembly_ftp_path(taxid, email, silent=True)
    except Exception:
        return
    if not ftp_path:
        return
    warn_gzip_size_vs_ncbi(path, identity=str(taxid), ftp_path=ftp_path)


def _raw_parent_for_dest(dest: Path, raw_dir: Optional[str] = None) -> Path:
    """Directory for NCBI ``*_genomic.fna.gz`` (sibling ``raw/``, not processed)."""
    if raw_dir:
        parent = Path(raw_dir)
        parent.mkdir(parents=True, exist_ok=True)
        return parent
    dest = Path(dest)
    if dest.name == "processed":
        parent = dest.parent / "raw"
    else:
        parent = dest / "raw"
    parent.mkdir(parents=True, exist_ok=True)
    return parent


def _ncbi_raw_cache_dir(output_path: Path) -> Tuple[Path, Path]:
    """Return ``(search_root, raw_dir)`` for taxid NCBI gzip downloads."""
    cache = genome_download_dir()
    try:
        cache.mkdir(parents=True, exist_ok=True)
        raw_dir = cache / "raw"
        raw_dir.mkdir(parents=True, exist_ok=True)
        return cache, raw_dir
    except OSError:
        raw_dir = Path(output_path) / "raw"
        try:
            raw_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            raw_dir = Path(output_path)
        return Path(output_path), raw_dir


def skip_new_ncbi_download(
    identity: str,
    *,
    skip_list: Optional[Collection[str]] = None,
    max_genome_mb: Any = None,
    ftp_path: Optional[str] = None,
    extra_ids: Sequence[str] = (),
    silent: bool = False,
) -> bool:
    """True when a *new* NCBI fetch should be skipped (skip-list or size cap).

    Already-on-disk genomes are the caller's job: this helper is only for the
    download branch. Default ``max_genome_mb`` is unlimited.
    """
    skip = (
        parse_genome_skip_list(skip_list, default_from_env=False)
        if skip_list is not None
        else parse_genome_skip_list(None, default_from_env=True)
    )
    max_mb = parse_max_genome_mb(max_genome_mb, default_from_env=max_genome_mb is None)
    ids = [identity, *extra_ids]
    for item in ids:
        if genome_id_in_skip_list(str(item), skip):
            if not silent:
                logger.info(
                    "Skipping NCBI download for %s (genome-skip-list)",
                    identity,
                )
            return True
    if max_mb == UNLIMITED_GENOME_MB or not ftp_path:
        return False
    size_mb = _remote_fasta_size_mb(ftp_path)
    if size_mb is not None and size_mb > max_mb:
        if not silent:
            logger.info(
                "Skipping NCBI download for %s: %.1f MB > max_genome_mb=%.1f",
                identity,
                size_mb,
                max_mb,
            )
        return True
    return False


def _effective_fetch_limits(
    max_genome_mb: Any = None,
    genome_skip_list: Any = None,
) -> Tuple[float, FrozenSet[str]]:
    return (
        parse_max_genome_mb(max_genome_mb, default_from_env=max_genome_mb is None),
        parse_genome_skip_list(genome_skip_list, default_from_env=genome_skip_list is None),
    )


def _accession_from_ftp_path(ftp_path: str) -> str:
    name = os.path.basename(str(ftp_path or "").rstrip("/"))
    parts = name.split("_")
    if len(parts) >= 2:
        cand = f"{parts[0]}_{parts[1]}"
        if is_assembly_accession(cand):
            return cand
    return ""


def _assembly_doc_for_taxid(taxid: str | int, email: str, silent: bool = False) -> Optional[dict]:
    """NCBI assembly DocumentSummary for ``taxid`` (descendant complete RefSeq)."""
    taxid = _normalize_taxid(taxid)
    Entrez.email = email
    search_term = f"txid{taxid}[Organism:exp] AND \"complete genome\"[filter] AND \"latest refseq\"[filter]"
    if not silent:
        logger.info(f"Searching with term: {search_term}")

    handle = Entrez.esearch(db="assembly", term=search_term, retmax=1)
    record = Entrez.read(handle)
    handle.close()
    if not record["IdList"]:
        if not silent:
            logger.warning(f"No genome found for taxid {taxid}")
        return None

    assembly_id = record["IdList"][0]
    handle = Entrez.esummary(db="assembly", id=assembly_id)
    summary = Entrez.read(handle)
    handle.close()
    docs = summary["DocumentSummarySet"]["DocumentSummary"]
    return docs[0] if docs else None


def _assembly_ftp_path(taxid: str | int, email: str, silent: bool = False) -> Optional[str]:
    """Return RefSeq (preferred) or GenBank FTP path for a taxid assembly."""
    doc = _assembly_doc_for_taxid(taxid, email, silent=silent)
    if not doc:
        return None
    ftp_path = doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank")
    if not ftp_path and not silent:
        logger.warning(f"No FTP path found for taxid {taxid}")
    return ftp_path or None


def _assembly_record(accession: str, email: str, silent: bool = False) -> Optional[dict]:
    accession = normalize_accession(accession) or accession
    Entrez.email = email
    handle = Entrez.esearch(db="assembly", term=f"{accession}[Assembly Accession]", retmax=1)
    record = Entrez.read(handle)
    handle.close()
    if not record.get("IdList"):
        if not silent:
            logger.warning("No NCBI assembly for %s", accession)
        return None
    handle = Entrez.esummary(db="assembly", id=record["IdList"][0])
    summary = Entrez.read(handle)
    handle.close()
    docs = summary["DocumentSummarySet"]["DocumentSummary"]
    return docs[0] if docs else None


def _assembly_ftp_for_accession(accession: str, email: str, silent: bool = False) -> Optional[str]:
    doc = _assembly_record(accession, email, silent=silent)
    if not doc:
        return None
    return doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank") or None


def fetch_assembly_processed(
    accession: str,
    dest_dir: str,
    email: str,
    *,
    silent: bool = False,
    keep_raw: bool = False,
    raw_dir: Optional[str] = None,
    index: bool = True,
    max_genome_mb: Any = None,
    genome_skip_list: Any = None,
) -> Optional[str]:
    """Download NCBI assembly ``accession`` and write ``{accession}.fa.gz``."""
    accession = normalize_accession(accession) or accession
    indexed = resolve_indexed_path(accession)
    if indexed is not None:
        if not silent:
            logger.info("Using indexed genome %s: %s", accession, indexed)
        return str(indexed)
    dest = Path(dest_dir)
    dest.mkdir(parents=True, exist_ok=True)
    existing = dest / processed_filename(accession)
    taxid, species = assembly_taxonomy(accession, email, silent=True)
    if existing.is_file():
        if index:
            index_processed_file(
                existing,
                accession=accession,
                move_to=dest,
                taxid=taxid,
                species_taxid=species,
            )
        return str(existing)
    ftp_path = _assembly_ftp_for_accession(accession, email, silent=silent)
    if not ftp_path:
        return None
    if skip_new_ncbi_download(
        accession,
        skip_list=genome_skip_list,
        max_genome_mb=max_genome_mb,
        ftp_path=ftp_path,
        extra_ids=[taxid, species],
        silent=silent,
    ):
        return None
    asm_name = os.path.basename(ftp_path)
    url = f"{ftp_path}/{asm_name}_genomic.fna.gz".replace("ftp://", "https://")
    raw_parent = _raw_parent_for_dest(dest, raw_dir)
    raw_path = raw_parent / raw_filename(accession)
    try:
        if raw_path.is_file():
            warn_gzip_size_vs_ncbi(raw_path, identity=accession, ftp_path=ftp_path)
            if not silent:
                logger.info("Reusing cached NCBI gzip %s", raw_path)
        else:
            if not silent:
                logger.info("Downloading genome from %s", url)
            _download_url(url, raw_path)
    except Exception as exc:
        if not silent:
            logger.error("Failed to download %s: %s", accession, exc)
        return None
    processed = dest / processed_filename(accession)
    try:
        preprocess_fasta(
            input_file=str(raw_path),
            output_file=str(processed),
            mutation_rate=0.0,
            include_percent=100.0,
        )
    except Exception as exc:
        logger.error("Failed to process %s: %s", accession, exc)
        return None
    if not keep_raw:
        try:
            raw_path.unlink()
        except OSError:
            pass
    if index:
        processed = index_processed_file(
            processed,
            accession=accession,
            move_to=dest,
            taxid=taxid,
            species_taxid=species,
        )
    return str(processed)


def assembly_taxonomy(accession: str, email: str = "", silent: bool = True) -> Tuple[str, str]:
    """Return ``(taxid, species_taxid)`` for an assembly accession."""
    try:
        doc = _assembly_record(
            accession, email or default_entrez_email(), silent=silent
        )
    except Exception:
        return "", ""
    if not doc:
        return "", ""
    taxid = str(doc.get("Taxid") or "").strip()
    species = str(doc.get("SpeciesTaxid") or taxid).strip()
    return taxid, species


def assembly_taxid(accession: str, email: str = "", silent: bool = True) -> str:
    """NCBI Taxid for an assembly accession (empty if lookup fails)."""
    taxid, _species = assembly_taxonomy(accession, email, silent=silent)
    return taxid


def materialize_accessions(
    accessions: Sequence[str],
    *,
    output_dir: str,
    email: str,
    reindex_mode: int = 1,
    keep_raw: bool = False,
    silent: bool = False,
    max_genome_mb: Any = None,
    genome_skip_list: Any = None,
) -> List[str]:
    """Download/process accessions according to ``samovar generate --reindex``.

    0: ``$out/.genomes/processed``, do not index
    1: ``samovar_database/processed`` and index
    2: ``$out/.genomes/processed`` and index in place
    """
    mode = int(reindex_mode)
    out = Path(output_dir)
    run_dest = run_processed_dir(out)
    if mode == 1:
        dest = processed_storage_dir()
        raw_dest = raw_storage_dir()
        do_index = True
    else:
        dest = run_dest
        raw_dest = run_raw_dir(out)
        do_index = mode == 2
    dest.mkdir(parents=True, exist_ok=True)
    run_dest.mkdir(parents=True, exist_ok=True)
    paths: List[str] = []
    for acc in accessions:
        acc = normalize_accession(acc) or str(acc).strip()
        if not acc:
            continue
        indexed = resolve_indexed_path(acc)
        if indexed is not None:
            if not silent:
                logger.info("Skip download, indexed: %s", indexed)
            paths.append(str(stage_into_dir(indexed, run_dest)))
            continue
        local = run_dest / processed_filename(acc)
        already = dest / processed_filename(acc)
        source = None
        if already.is_file():
            source = already
        elif local.is_file():
            source = local
        if source is not None:
            if not silent:
                logger.info("Reuse existing processed genome %s: %s", acc, source)
            if do_index:
                source = index_processed_file(source, accession=acc, move_to=dest)
            paths.append(str(stage_into_dir(source, run_dest)))
            continue
        fetched = fetch_assembly_processed(
            acc,
            str(dest),
            email,
            silent=silent,
            keep_raw=keep_raw,
            raw_dir=str(raw_dest),
            index=do_index,
            max_genome_mb=max_genome_mb,
            genome_skip_list=genome_skip_list,
        )
        if fetched:
            paths.append(str(stage_into_dir(fetched, run_dest)))
        else:
            logger.warning("Could not materialize %s", acc)
    return paths


def _library_annotation_candidates(taxid: str, dest_name: str) -> List[Path]:
    taxid = _normalize_taxid(taxid)
    roots: List[Path] = []
    try:
        roots.append(proteome_storage_dir())
    except Exception:
        pass
    indexed = resolve_indexed_path(taxid)
    if indexed is not None:
        roots.append(indexed.parent)
        if indexed.parent.name == "processed":
            roots.append(indexed.parent.parent / "proteomes")
    names = {dest_name, f"{taxid}{Path(dest_name).suffix}"}
    if dest_name.endswith(".gz"):
        names.add(dest_name[:-3])
    found: List[Path] = []
    for root in roots:
        for name in names:
            cand = root / name
            if cand.is_file():
                found.append(cand)
    return found


def _download_assembly_file(
    taxid: str | int,
    output_folder: str,
    email: str,
    suffix: str,
    dest_name: str,
    silent: bool = False,
) -> Optional[str]:
    """Download an NCBI assembly file (suffix like '_protein.faa.gz').

    Files stay gzipped on disk; callers that need plain text open them via seqio.
    Shared copies live under ``samovar_database/proteomes``.
    """
    taxid = _normalize_taxid(taxid)
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    dest = output_path / dest_name
    dest_gz = dest if is_gzip_path(dest) else Path(str(dest) + ".gz")
    dest_plain = Path(str(dest).removesuffix(".gz")) if is_gzip_path(dest) else dest
    for candidate in (dest, dest_gz, dest_plain, *_library_annotation_candidates(taxid, dest_name)):
        if candidate.exists():
            return str(candidate)
    try:
        lib = proteome_storage_dir()
        lib.mkdir(parents=True, exist_ok=True)
        lib_dest = lib / dest_name
        if lib_dest.exists():
            return str(lib_dest)
    except Exception:
        lib_dest = dest_gz
    try:
        ftp_path = _assembly_ftp_path(taxid, email, silent=silent)
        if not ftp_path:
            return None
        asm_name = os.path.basename(ftp_path)
        url = f"{ftp_path}/{asm_name}{suffix}".replace("ftp://", "https://")
        if not silent:
            logger.info(f"Downloading {url}")
        _download_url(url, lib_dest if lib_dest.parent.exists() else dest_gz)
        return str(lib_dest if lib_dest.exists() else dest_gz)
    except Exception as e:
        if not silent:
            logger.error(f"Failed to download {suffix} for taxid {taxid}: {e}")
        return None


def fetch_proteome(
    taxid: str | int,
    output_folder: str,
    email: str,
    silent: bool = False,
) -> Optional[str]:
    """Download NCBI protein FASTA (`_protein.faa.gz`) for a taxid."""
    return _download_assembly_file(
        taxid,
        output_folder,
        email,
        "_protein.faa.gz",
        f"{_normalize_taxid(taxid)}.faa.gz",
        silent,
    )


def fetch_gff(
    taxid: str | int,
    output_folder: str,
    email: str,
    silent: bool = False,
) -> Optional[str]:
    """Download NCBI genomic GFF (`_genomic.gff.gz`) for a taxid."""
    return _download_assembly_file(
        taxid,
        output_folder,
        email,
        "_genomic.gff.gz",
        f"{_normalize_taxid(taxid)}.gff.gz",
        silent,
    )


def fetch_genome_raw(
    taxid: str|int,
    output_folder: str,
    email: str,
    reference_only: bool = True,
    silent: bool = False,
    reuse: Optional[bool] = None,
    max_genome_mb: Any = None,
    genome_skip_list: Any = None,
    ) -> Optional[str]:
    """Fetch a genome for ``taxid``: run dir, NCBI library cache, then NCBI.

    Bundled ``data/test_genomes`` stubs are not used unless explicitly allowed.
    New downloads land in the install ``genomes/raw`` cache; processed FASTA
    is written later by ``fetch_genome``.
    """
    taxid = _normalize_taxid(taxid)
    Entrez.email = email

    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    possible_exts = [".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"]
    for ext in possible_exts:
        candidate = output_path / f"{taxid}{ext}"
        if candidate.exists():
            _warn_existing_ncbi_gzip(candidate, taxid, email)
            if not silent:
                logger.info(f"Genome for taxid {taxid} already exists at {candidate}")
            return str(candidate)

    do_reuse = reuse_enabled(reuse)
    if do_reuse:
        library = find_library_genome(taxid, extra=[output_path], include_test_genomes=False)
        if library is not None:
            placed = place_genome(library, output_path, taxid)
            if not silent:
                logger.info(f"Reusing library genome for taxid {taxid}: {placed}")
            return str(placed)

    if allow_test_genomes():
        bundled = _copy_bundled_genome(taxid, output_folder)
        if bundled is not None:
            if not silent:
                logger.warning(
                    "Using truncated bundled test genome for taxid %s: %s "
                    "(SAMOVAR_ALLOW_TEST_GENOMES / --test-genomes)",
                    taxid,
                    bundled,
                )
            return bundled

    cache, raw_dir = _ncbi_raw_cache_dir(output_path)
    genome_path = raw_dir / f"{taxid}.fna.gz"
    existing_cache = None
    for folder in (raw_dir, cache):
        for ext in possible_exts:
            candidate = folder / f"{taxid}{ext}"
            if candidate.exists():
                existing_cache = candidate
                break
        if existing_cache is not None:
            break
    if existing_cache is not None:
        genome_path = existing_cache
        _warn_existing_ncbi_gzip(genome_path, taxid, email)
    else:
        old_timeout = socket.getdefaulttimeout()
        socket.setdefaulttimeout(45)
        try:
            max_mb, skip = _effective_fetch_limits(max_genome_mb, genome_skip_list)
            extra_ids: List[str] = []
            if genome_id_in_skip_list(taxid, skip):
                if not silent:
                    logger.info("Skipping NCBI download for %s (genome-skip-list)", taxid)
                return None
            ftp_path = _assembly_ftp_path(taxid, email, silent=silent)
            if not ftp_path:
                return None
            try:
                doc = _assembly_doc_for_taxid(taxid, email, silent=True)
                if doc:
                    extra_ids.append(str(doc.get("AssemblyAccession") or ""))
                    extra_ids.append(str(doc.get("Taxid") or ""))
            except Exception:
                pass
            if skip_new_ncbi_download(
                taxid,
                skip_list=skip,
                max_genome_mb=max_mb,
                ftp_path=ftp_path,
                extra_ids=extra_ids,
                silent=silent,
            ):
                return None
            http_url = _genomic_fna_gz_url(ftp_path)
            if not silent:
                logger.info(f"Downloading genome from {http_url}")
            _download_url(http_url, genome_path)
        except Exception as e:
            if not silent:
                logger.error(f"Error fetching genome for taxid {taxid}: {str(e)}")
            return None
        finally:
            socket.setdefaulttimeout(old_timeout)

    ensure_genome_config()
    register_genome_dir(genome_path.parent)
    return str(genome_path)

def fetch_genome(
    taxid: str|int,
    output_folder: str,
    email: str,
    reference_only: bool = True,
    silent: bool = False,
    gzip_genomes: bool = True,
    reuse: Optional[bool] = None,
    index: bool = True,
    max_genome_mb: Any = None,
    genome_skip_list: Any = None,
    ) -> Optional[str]:
    """Fetch/process a genome. Size and skip-list apply only to new NCBI downloads.

    Local, indexed, and library copies are reused at any size (including skip-list
    taxids that are already on disk).
    """
    max_mb, skip = _effective_fetch_limits(max_genome_mb, genome_skip_list)

    if is_assembly_accession(str(taxid)):
        return fetch_assembly_processed(
            str(taxid),
            output_folder,
            email,
            silent=silent,
            keep_raw=False,
            index=index,
            max_genome_mb=max_mb,
            genome_skip_list=skip,
        )

    if isinstance(taxid, str) and taxid.split(".")[0].isdigit():
        taxid = taxid.split(".")[0]

    existing = find_existing_processed_genome(output_folder, taxid)
    if existing is not None:
        if gzip_genomes and not is_gzip_path(existing):
            existing = gzip_file(existing)
        elif not gzip_genomes and is_gzip_path(existing):
            existing = gunzip_file(existing, remove_source=True)
        if not silent:
            logger.info(f"Processed genome for taxid {taxid} already exists at {existing}")
        return str(existing)

    indexed = resolve_indexed_path(str(taxid))
    if indexed is not None:
        if not silent:
            logger.info("Using indexed genome for %s: %s", taxid, indexed)
        return str(indexed)

    do_reuse = reuse_enabled(reuse)
    if do_reuse:
        library = find_library_genome(
            taxid, extra=[output_folder], include_test_genomes=False, processed_only=True
        )
        if library is not None:
            placed = place_genome(library, output_folder, taxid)
            if gzip_genomes and not is_gzip_path(placed):
                placed = Path(gzip_file(placed))
            elif not gzip_genomes and is_gzip_path(placed):
                placed = Path(gunzip_file(placed, remove_source=False))
                placed = place_genome(placed, output_folder, taxid)
            if index:
                try:
                    catalog_processed_genome(
                        placed,
                        taxid=str(taxid),
                        keep_src=True,
                        stage_dir=output_folder,
                    )
                except Exception:
                    logger.warning("Could not catalog reused genome for %s", taxid, exc_info=True)
            if not silent:
                logger.info(f"Reusing processed library genome for taxid {taxid}: {placed}")
            return str(placed)

    genome_path = fetch_genome_raw(
        taxid,
        output_folder,
        email,
        reference_only,
        silent,
        reuse=reuse,
        max_genome_mb=max_mb,
        genome_skip_list=skip,
    )
    if genome_path is None:
        return None

    processed_dir = cache_processed_dir()
    try:
        processed_dir.mkdir(parents=True, exist_ok=True)
    except OSError:
        # Never fall back to $HOME — keep artifacts next to the run genomes.
        processed_dir = Path(output_folder)
        try:
            processed_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            pass
    cache_processed = processed_genome_path(processed_dir, taxid, gzip_genomes=gzip_genomes)
    run_processed = processed_genome_path(output_folder, taxid, gzip_genomes=gzip_genomes)

    try:
        preprocess_fasta(
            input_file=genome_path,
            output_file=str(cache_processed),
            mutation_rate=0.0,
            include_percent=100.0,
        )
        ensure_genome_config()
        register_genome_dir(processed_dir)
        if Path(cache_processed).resolve() != Path(run_processed).resolve():
            place_genome(cache_processed, output_folder, taxid)
        logger.info(f"Successfully processed genome for taxid {taxid}")
        found = find_existing_processed_genome(output_folder, taxid)
        result = str(found or cache_processed)
        if index:
            acc = accession_from_fasta(genome_path) or accession_from_fasta(result)
            species = ""
            try:
                _tax, species = assembly_taxonomy(acc, email, silent=True) if is_assembly_accession(acc) else ("", "")
            except Exception:
                species = ""
            if not acc or acc == str(taxid):
                try:
                    doc = _assembly_doc_for_taxid(taxid, email, silent=True)
                    if doc:
                        acc = str(doc.get("AssemblyAccession") or acc or taxid)
                        species = str(doc.get("SpeciesTaxid") or species or taxid)
                except Exception:
                    pass
            try:
                catalog_processed_genome(
                    result,
                    taxid=str(taxid),
                    accession=acc if is_assembly_accession(acc) else "",
                    species_taxid=species or str(taxid),
                    keep_src=True,
                    stage_dir=output_folder,
                )
            except Exception:
                logger.warning("Processed %s but could not update install config", taxid, exc_info=True)
        return result
    except Exception as e:
        logger.error(f"Failed to process genome for taxid {taxid}: {str(e)}")
        return None


def generate_random_taxids(
    group: str = "Bacteria",
    N: int = 10,
    silent: bool = False,
    max_genome_mb: float | None = None,
    genome_skip_list: Any = None,
    oversample: int = 5,
) -> list[str]:
    """
    Generate a list of random taxids for a given taxonomic group.
    
    Args:
        group (str): Taxonomic group to sample from (default: "Bacteria")
        N (int): Number of unique taxids to generate
        silent (bool): If True, suppress logging output
        max_genome_mb: Prefer assemblies whose genomic FASTA is under this size
            (via HTTP HEAD Content-Length). Default unlimited (``inf``).
            ``None`` reads ``SAMOVAR_MAX_GENOME_MB``.
        genome_skip_list: Taxids/accessions to ignore when sampling.
        oversample: Collect up to N * oversample RefSeq candidates before sampling
        
    Returns:
        list[str]: List of unique taxids
    """
    # Set up Entrez
    if not hasattr(Entrez, 'email'):
        raise ValueError("Entrez.email must be set before calling this function")

    max_mb, skip = _effective_fetch_limits(max_genome_mb, genome_skip_list)
    
    old_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(45)
    try:
        # Search for organisms with complete genomes
        search_term = f'"{group}"[Organism] AND "latest refseq"[filter] AND "complete genome"[filter]'
        if not silent:
            logger.info(f"Searching with term: {search_term}")
        
        def search_func():
            handle = Entrez.esearch(db="assembly", term=search_term, retmax=1000)
            record = Entrez.read(handle)
            handle.close()
            return record
            
        record = _entrez_retry(search_func)
        
        if not silent:
            logger.info(f"Found {len(record['IdList'])} assemblies")
        
        if not record["IdList"]:
            if not silent:
                logger.warning(f"No genomes found for group {group}")
            return []
            
        # Shuffle the assembly IDs to get more random results
        random.shuffle(record["IdList"])
            
        # Collect (size_mb, taxid) candidates; prefer smaller genomes when filtered
        candidates: list[tuple[float, str]] = []
        target = max(N * oversample, N)
        for assembly_id in record["IdList"]:
            def summary_func():
                handle = Entrez.esummary(db="assembly", id=assembly_id)
                summary = Entrez.read(handle)
                handle.close()
                return summary
                
            summary = _entrez_retry(summary_func)
            
            # Check if this assembly has a RefSeq FTP path
            doc_summary = summary["DocumentSummarySet"]["DocumentSummary"][0]
            ftp_path = doc_summary.get("FtpPath_RefSeq")
            if not ftp_path:
                continue
            taxid = doc_summary.get("Taxid")
            if not taxid:
                continue
            taxid = str(taxid)
            if genome_id_in_skip_list(taxid, skip):
                if not silent:
                    logger.info(f"Skipping taxid {taxid} (genome-skip-list)")
                continue
            size_mb = float("inf")
            if max_mb != UNLIMITED_GENOME_MB:
                size = _remote_fasta_size_mb(ftp_path)
                if size is None:
                    size_mb = float("inf")
                else:
                    size_mb = size
                if size_mb != float("inf") and size_mb > max_mb:
                    if not silent:
                        logger.info(f"Skipping taxid {taxid}: {size_mb:.1f} MB > {max_mb} MB")
                    continue
            if any(t == taxid for _, t in candidates):
                continue
            candidates.append((size_mb, taxid))
            if not silent:
                logger.info(f"Found taxid: {taxid}" + (f" ({size_mb:.1f} MB)" if size_mb != float("inf") else ""))
            if len(candidates) >= target:
                break
        
        if not candidates:
            if not silent:
                logger.warning(f"No suitable genomes found for group {group}")
            return []

        candidates.sort(key=lambda x: x[0])
        taxids = [t for _, t in candidates[: max(N * 2, N)]]
        if len(taxids) > N:
            # Prefer smaller half, then sample
            taxids = taxids[:N]
            
        if not silent:
            logger.info(f"Returning {len(taxids)} taxids")
        return taxids
        
    except Exception as e:
        if not silent:
            logger.error(f"Error generating random taxids: {str(e)}")
        return []
    finally:
        socket.setdefaulttimeout(old_timeout)


def main():
    """Main function to process genomes from random taxids."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Process genomes from random taxids')
    parser.add_argument('--group', type=str, default='Bacteria',
                      help='Taxonomic group to sample from (default: Bacteria)')
    parser.add_argument('--N', type=int, default=10,
                      help='Number of genomes to process (default: 10)')
    parser.add_argument('--email', type=str, default=None,
                      help='Email for NCBI Entrez (default: NCBI_EMAIL / ENTREZ_EMAIL / SAMOVAR_EMAIL, else anonymous@example.com)')
    parser.add_argument('--output-dir', type=str, default='genomes',
                      help='Output directory for genomes (default: genomes)')
    parser.add_argument(
        '--accessions',
        nargs='+',
        default=None,
        help='Download these NCBI assembly accessions as {accession}.fa.gz',
    )
    parser.add_argument(
        '--reindex',
        type=int,
        default=0,
        choices=(0, 1, 2),
        help='With --accessions: 0 run .genomes only; 1 samovar_database+index; 2 run+index',
    )
    parser.add_argument('--silent', action='store_true',
                      help='Suppress logging output and show progress bars instead')
    parser.add_argument(
        '--max-genome-mb',
        type=argparse_max_genome_mb,
        default=UNLIMITED_GENOME_MB,
        help='Skip new NCBI assemblies larger than this many MB (default: inf; 0 also unlimited)',
    )
    parser.add_argument(
        '--genome-skip-list',
        dest='genome_skip_list',
        default=None,
        help='Comma-separated taxids/accessions to skip on new NCBI downloads',
    )
    parser.add_argument(
        '--gzip-genomes',
        dest='gzip_genomes',
        action='store_true',
        default=True,
        help='Keep processed genomes gzip-compressed (default)',
    )
    parser.add_argument(
        '--no-gzip-genomes',
        dest='gzip_genomes',
        action='store_false',
        help='Write uncompressed processed FASTA',
    )
    
    args = parser.parse_args()
    args.email = args.email or default_entrez_email()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    Entrez.email = args.email

    if args.accessions:
        materialize_accessions(
            args.accessions,
            output_dir=str(output_dir),
            email=args.email,
            reindex_mode=int(args.reindex),
            keep_raw=False,
            silent=args.silent,
            max_genome_mb=args.max_genome_mb,
            genome_skip_list=args.genome_skip_list,
        )
        return

    max_mb = parse_max_genome_mb(args.max_genome_mb, default_from_env=False)
    
    # Step 1: Generate candidate taxids (oversample so failed downloads can be replaced)
    if not args.silent:
        logger.info(f"Generating {args.N} random taxids for group {args.group}")
    taxids = generate_random_taxids(
        group=args.group,
        N=max(args.N * 3, args.N),
        silent=args.silent,
        max_genome_mb=max_mb,
        genome_skip_list=args.genome_skip_list,
    )
    
    if not taxids:
        if not args.silent:
            logger.error("No taxids found. Exiting.")
        return
    
    # Step 2: Fetch genomes until N successes
    successes = 0
    for taxid in tqdm(taxids, desc="Processing genomes", disable=not args.silent):
        if successes >= args.N:
            break
        if not args.silent:
            logger.info(f"Processing taxid {taxid}")
        
        genome_path = fetch_genome(
            taxid,
            str(output_dir),
            args.email,
            silent=args.silent,
            gzip_genomes=args.gzip_genomes,
            max_genome_mb=max_mb,
            genome_skip_list=args.genome_skip_list,
        )
        if not genome_path:
            if not args.silent:
                logger.warning(f"Failed to fetch genome for taxid {taxid}")
            continue
        processed = find_existing_processed_genome(output_dir, taxid)
        if processed is not None:
            successes += 1
    
    # Step 3: Cleanup - remove all files except processed FASTA (plain or gzipped)
    if not args.silent:
        logger.info("Cleaning up intermediate files...")
    for file in output_dir.glob("*"):
        name = file.name
        keep = (
            name.endswith(".fa.gz")
            or name.endswith(".fa")
            or "-processed.fasta" in name
            or "-processed.fna" in name
            or "-processed.fa" in name
        )
        if keep:
            continue
        try:
            file.unlink()
            if not args.silent:
                logger.info(f"Removed {file}")
        except Exception as e:
            if not args.silent:
                logger.error(f"Failed to remove {file}: {str(e)}")
    
    if not args.silent:
        logger.info(f"Processing complete! Got {successes}/{args.N} genomes.")

if __name__ == '__main__':
    main()
