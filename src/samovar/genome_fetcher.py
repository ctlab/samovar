"""
Genome fetching and taxonomy parsing functionality
"""

import os
import logging
import shutil
import socket
from typing import Optional
import urllib.request
from Bio import Entrez
from pathlib import Path
import random
from tqdm import tqdm
import time
from samovar.fasta_processor import preprocess_fasta
from samovar.seqio import (
    find_existing_processed_genome,
    gunzip_file,
    gzip_file,
    is_gzip_path,
    processed_genome_path,
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

    Used before NCBI so ISS / tests / air-gapped installs do not download
    reference assemblies (especially taxid 9606).
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
    with urllib.request.urlopen(req, timeout=timeout) as resp, open(dest, "wb") as out:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)


def _normalize_taxid(taxid: str | int) -> str:
    taxid = str(taxid)
    return taxid.split(".")[0]


def _assembly_ftp_path(taxid: str | int, email: str, silent: bool = False) -> Optional[str]:
    """Return RefSeq (preferred) or GenBank FTP path for a taxid assembly."""
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
    doc = summary["DocumentSummarySet"]["DocumentSummary"][0]
    ftp_path = doc.get("FtpPath_RefSeq") or doc.get("FtpPath_GenBank")
    if not ftp_path and not silent:
        logger.warning(f"No FTP path found for taxid {taxid}")
    return ftp_path or None


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
    """
    taxid = _normalize_taxid(taxid)
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    dest = output_path / dest_name
    dest_gz = dest if is_gzip_path(dest) else Path(str(dest) + ".gz")
    dest_plain = Path(str(dest).removesuffix(".gz")) if is_gzip_path(dest) else dest
    for candidate in (dest, dest_gz, dest_plain):
        if candidate.exists():
            return str(candidate)
    try:
        ftp_path = _assembly_ftp_path(taxid, email, silent=silent)
        if not ftp_path:
            return None
        asm_name = os.path.basename(ftp_path)
        url = f"{ftp_path}/{asm_name}{suffix}".replace("ftp://", "https://")
        if not silent:
            logger.info(f"Downloading {url}")
        _download_url(url, dest_gz)
        return str(dest_gz)
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
    silent: bool = False
    ) -> Optional[str]:
    """
    Fetch genome from NCBI for a given taxid.
    """
    taxid = _normalize_taxid(taxid)
    Entrez.email = email

    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)

    possible_exts = [".fa.gz", ".fna.gz", ".fasta.gz", ".fa", ".fna", ".fasta"]
    for ext in possible_exts:
        candidate = output_path / f"{taxid}{ext}"
        if candidate.exists():
            if not silent:
                logger.info(f"Genome for taxid {taxid} already exists at {candidate}")
            return str(candidate)

    bundled = _copy_bundled_genome(taxid, output_folder)
    if bundled is not None:
        if not silent:
            logger.info(f"Using bundled test genome for taxid {taxid}: {bundled}")
        return bundled

    genome_path = output_path / f"{taxid}.fna.gz"
    old_timeout = socket.getdefaulttimeout()
    socket.setdefaulttimeout(45)
    try:
        ftp_path = _assembly_ftp_path(taxid, email, silent=silent)
        if not ftp_path:
            return None
        asm_name = os.path.basename(ftp_path)
        http_url = f"{ftp_path}/{asm_name}_genomic.fna.gz".replace("ftp://", "https://")
        if not silent:
            logger.info(f"Downloading genome from {http_url}")
        _download_url(http_url, genome_path)
        return str(genome_path)
    except Exception as e:
        if not silent:
            logger.error(f"Error fetching genome for taxid {taxid}: {str(e)}")
        return None
    finally:
        socket.setdefaulttimeout(old_timeout)

def fetch_genome(
    taxid: str|int,
    output_folder: str,
    email: str,
    reference_only: bool = True,
    silent: bool = False,
    gzip_genomes: bool = True,
    ) -> Optional[str]:
    
    if isinstance(taxid, str):
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

    output_path = processed_genome_path(output_folder, taxid, gzip_genomes=gzip_genomes)

    genome_path = fetch_genome_raw(taxid, output_folder, email, reference_only, silent)
    if genome_path is None:
        return None

    try:
        preprocess_fasta(
            input_file=genome_path,
            output_file=str(output_path),
            mutation_rate=0.0,
            include_percent=100.0,
        )
        logger.info(f"Successfully processed genome for taxid {taxid}")
        return str(output_path)
    except Exception as e:
        logger.error(f"Failed to process genome for taxid {taxid}: {str(e)}")
        return None

def generate_random_taxids(
    group: str = "Bacteria",
    N: int = 10,
    silent: bool = False,
    max_genome_mb: float | None = 100.0,
    oversample: int = 5,
) -> list[str]:
    """
    Generate a list of random taxids for a given taxonomic group.
    
    Args:
        group (str): Taxonomic group to sample from (default: "Bacteria")
        N (int): Number of unique taxids to generate
        silent (bool): If True, suppress logging output
        max_genome_mb: Prefer assemblies whose genomic FASTA is under this size
            (via HTTP HEAD Content-Length). None disables size filtering.
        oversample: Collect up to N * oversample RefSeq candidates before sampling
        
    Returns:
        list[str]: List of unique taxids
    """
    # Set up Entrez
    if not hasattr(Entrez, 'email'):
        raise ValueError("Entrez.email must be set before calling this function")
    
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
            size_mb = float("inf")
            if max_genome_mb is not None:
                try:
                    asm_name = os.path.basename(ftp_path)
                    url = f"{ftp_path}/{asm_name}_genomic.fna.gz".replace("ftp://", "https://")
                    req = urllib.request.Request(url, method="HEAD")
                    with urllib.request.urlopen(req, timeout=30) as resp:
                        cl = resp.headers.get("Content-Length")
                        if cl:
                            size_mb = int(cl) / (1024 * 1024)
                except Exception:
                    # If size unknown, keep candidate but rank last
                    size_mb = float("inf")
                if size_mb != float("inf") and size_mb > max_genome_mb:
                    if not silent:
                        logger.info(f"Skipping taxid {taxid}: {size_mb:.1f} MB > {max_genome_mb} MB")
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
    parser.add_argument('--silent', action='store_true',
                      help='Suppress logging output and show progress bars instead')
    parser.add_argument('--max-genome-mb', type=float, default=100.0,
                      help='Skip assemblies larger than this many MB (0 disables)')
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
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up Entrez email
    Entrez.email = args.email
    max_mb = None if args.max_genome_mb <= 0 else args.max_genome_mb
    
    # Step 1: Generate candidate taxids (oversample so failed downloads can be replaced)
    if not args.silent:
        logger.info(f"Generating {args.N} random taxids for group {args.group}")
    taxids = generate_random_taxids(
        group=args.group,
        N=max(args.N * 3, args.N),
        silent=args.silent,
        max_genome_mb=max_mb,
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
        keep = "-processed.fasta" in name or "-processed.fna" in name or "-processed.fa" in name
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
