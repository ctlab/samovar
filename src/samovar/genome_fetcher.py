"""
Genome fetching and taxonomy parsing functionality
"""

import os
import logging
from typing import Optional
import urllib.request
from Bio import Entrez
from pathlib import Path
import random
import gzip
import shutil
from tqdm import tqdm
import time

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

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
            if "429" in str(e) or "Too Many Requests" in str(e):
                if attempt < max_retries - 1:
                    logger.warning(f"Rate limited, retrying in {delay} seconds...")
                    time.sleep(delay)
                    delay *= 2  # Exponential backoff
                continue
            raise
    
    raise last_exception

def fetch_genome(
    taxid: str|int,
    output_folder: str,
    email: str,
    reference_only: bool = True,
    silent: bool = False
    ) -> Optional[str]:
    """
    Fetch genome from NCBI for a given taxid.
    
    Args:
        taxid (str): NCBI taxonomy ID
        output_folder (str): Path to output folder
        email (str): Email for NCBI Entrez
        reference_only (bool): If True, only fetch reference genome
        silent (bool): If True, suppress logging output
        
    Returns:
        Optional[str]: Path to downloaded genome file or None if failed
    """
    if isinstance(taxid, str):
        taxid = taxid.split(".")[0]
    # Set up Entrez
    Entrez.email = email
    
    # Create output folder if it doesn't exist
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Check if genome already exists
    # Try to match any of .fa.gz, .fna.gz, .fasta.gz
    possible_exts = [".fa.gz", ".fna.gz", ".fasta.gz"]
    for ext in possible_exts:
        candidate = output_path / f"{taxid}{ext}"
        if candidate.exists():
            genome_path = candidate
            break
    else:
        # Default to .fna.gz for download if not present
        genome_path = output_path / f"{taxid}.fna.gz"
    if genome_path.exists():
        if not silent:
            logger.info(f"Genome for taxid {taxid} already exists at {genome_path}")
        return str(genome_path)
    
    try:
        # Search for genome assembly with complete genome
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
            
        # Get assembly details
        assembly_id = record["IdList"][0]
        handle = Entrez.esummary(db="assembly", id=assembly_id)
        summary = Entrez.read(handle)
        handle.close()
        
        # Get FTP path - prefer RefSeq path
        ftp_path = summary["DocumentSummarySet"]["DocumentSummary"][0].get("FtpPath_RefSeq")
        if not ftp_path:
            ftp_path = summary["DocumentSummarySet"]["DocumentSummary"][0].get("FtpPath_GenBank")
            
        if not ftp_path:
            if not silent:
                logger.warning(f"No FTP path found for taxid {taxid}")
            return None
            
        # Construct download URL
        asm_name = os.path.basename(ftp_path)
        url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
        
        # Convert FTP URL to HTTP URL
        http_url = url.replace('ftp://', 'https://')
        
        # Download and save
        try:
            if not silent:
                logger.info(f"Downloading genome from {http_url}")
            urllib.request.urlretrieve(http_url, genome_path)
            if not silent:
                logger.info(f"Successfully downloaded genome for taxid {taxid}")
            return str(genome_path)
        except Exception as e:
            if not silent:
                logger.error(f"Failed to download genome for taxid {taxid}: {str(e)}")
            return None
            
    except Exception as e:
        if not silent:
            logger.error(f"Error fetching genome for taxid {taxid}: {str(e)}")
        return None

def generate_random_taxids(group: str = "Bacteria", N: int = 10, silent: bool = False) -> list[str]:
    """
    Generate a list of random taxids for a given taxonomic group.
    
    Args:
        group (str): Taxonomic group to sample from (default: "Bacteria")
        N (int): Number of unique taxids to generate
        silent (bool): If True, suppress logging output
        
    Returns:
        list[str]: List of unique taxids
    """
    # Set up Entrez
    if not hasattr(Entrez, 'email'):
        raise ValueError("Entrez.email must be set before calling this function")
    
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
            
        # Get assembly details for each ID to extract taxids
        taxids = set()
        for assembly_id in record["IdList"]:
            def summary_func():
                handle = Entrez.esummary(db="assembly", id=assembly_id)
                summary = Entrez.read(handle)
                handle.close()
                return summary
                
            summary = _entrez_retry(summary_func)
            
            # Check if this assembly has a RefSeq FTP path
            doc_summary = summary["DocumentSummarySet"]["DocumentSummary"][0]
            if doc_summary.get("FtpPath_RefSeq"):
                taxid = doc_summary.get("Taxid")
                if taxid:
                    taxids.add(str(taxid))
                    if not silent:
                        logger.info(f"Found taxid: {taxid}")
            
            if len(taxids) >= N:
                break
        
        # Convert to list and take random sample
        taxids = list(taxids)
        if len(taxids) > N:
            taxids = random.sample(taxids, N)
            
        if not silent:
            logger.info(f"Returning {len(taxids)} taxids")
        return taxids
        
    except Exception as e:
        if not silent:
            logger.error(f"Error generating random taxids: {str(e)}")
        return []

def main():
    """Main function to process genomes from random taxids."""
    import argparse
    import subprocess
    from samovar.fasta_processor import preprocess_fasta
    
    parser = argparse.ArgumentParser(description='Process genomes from random taxids')
    parser.add_argument('--group', type=str, default='Bacteria',
                      help='Taxonomic group to sample from (default: Bacteria)')
    parser.add_argument('--N', type=int, default=10,
                      help='Number of genomes to process (default: 10)')
    parser.add_argument('--email', type=str, default="dsmutin@gmail.com",
                      help='Email for NCBI Entrez')
    parser.add_argument('--output-dir', type=str, default='genomes',
                      help='Output directory for genomes (default: genomes)')
    parser.add_argument('--silent', action='store_true',
                      help='Suppress logging output and show progress bars instead')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Set up Entrez email
    Entrez.email = args.email
    
    # Step 1: Generate random taxids
    if not args.silent:
        logger.info(f"Generating {args.N} random taxids for group {args.group}")
    taxids = generate_random_taxids(group=args.group, N=args.N, silent=args.silent)
    
    if not taxids:
        if not args.silent:
            logger.error("No taxids found. Exiting.")
        return
    
    # Step 2: Fetch genomes
    for taxid in tqdm(taxids, desc="Processing genomes", disable=not args.silent):
        if not args.silent:
            logger.info(f"Processing taxid {taxid}")
        
        # Fetch genome
        genome_path = fetch_genome(taxid, str(output_dir), args.email, silent=args.silent)
        if not genome_path:
            if not args.silent:
                logger.warning(f"Failed to fetch genome for taxid {taxid}")
            continue
            
        # Step 3: Un-gzip if needed
        if genome_path.endswith('.gz'):
            try:
                with gzip.open(genome_path, 'rb') as f_in, open(genome_path[:-3], 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                if not args.silent:
                    logger.info(f"Unzipped {genome_path} to {genome_path[:-3]}")
                genome_path = genome_path[:-3]
            except Exception as e:
                if not args.silent:
                    logger.error(f"Failed to unzip {genome_path}: {str(e)}")
                continue
        
        # Step 4: Preprocess FASTA
        output_path = output_dir / f"{taxid}-processed.fasta"
        try:
            preprocess_fasta(
                input_file=genome_path,
                output_file=str(output_path),
                mutation_rate=0.0,
                include_percent=100.0
            )
            if not args.silent:
                logger.info(f"Successfully processed genome for taxid {taxid}")
        except Exception as e:
            if not args.silent:
                logger.error(f"Failed to process genome for taxid {taxid}: {str(e)}")
    
    # Step 5: Cleanup - remove all files except processed FASTA files
    if not args.silent:
        logger.info("Cleaning up intermediate files...")
    for file in output_dir.glob("*"):
        if not file.name.endswith("-processed.fasta"):
            try:
                file.unlink()
                if not args.silent:
                    logger.info(f"Removed {file}")
            except Exception as e:
                if not args.silent:
                    logger.error(f"Failed to remove {file}: {str(e)}")
    
    if not args.silent:
        logger.info("Processing complete!")

if __name__ == '__main__':
    main()