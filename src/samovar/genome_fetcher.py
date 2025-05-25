"""
Genome fetching and taxonomy parsing functionality
"""

import os
import logging
from typing import Optional, List, Set
import urllib.request
import pandas as pd
from Bio import Entrez
from Bio import SeqIO
from pathlib import Path

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def fetch_genome(taxid: str, output_folder: str, email: str, reference_only: bool = True) -> Optional[str]:
    """
    Fetch genome from NCBI for a given taxid.
    
    Args:
        taxid (str): NCBI taxonomy ID
        output_folder (str): Path to output folder
        email (str): Email for NCBI Entrez
        reference_only (bool): If True, only fetch reference genome
        
    Returns:
        Optional[str]: Path to downloaded genome file or None if failed
    """
    # Set up Entrez
    Entrez.email = email
    
    # Create output folder if it doesn't exist
    output_path = Path(output_folder)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Check if genome already exists
    # Try to match any of .fa, .fna, .fasta
    possible_exts = [".fa", ".fna", ".fasta"]
    for ext in possible_exts:
        candidate = output_path / f"{taxid}{ext}"
        if candidate.exists():
            genome_path = candidate
            break
    else:
        # Default to .fna for download if not present
        genome_path = output_path / f"{taxid}.fna"
    if genome_path.exists():
        logger.info(f"Genome for taxid {taxid} already exists at {genome_path}")
        return str(genome_path)
    
    try:
        # Search for genome assembly - modified search term to be more inclusive
        search_term = f"txid{taxid}[Organism:exp]"
        if reference_only:
            search_term += " AND refseq[filter]"
            
        handle = Entrez.esearch(db="assembly", term=search_term, retmax=1)
        record = Entrez.read(handle)
        handle.close()
        
        if not record["IdList"]:
            logger.warning(f"No genome found for taxid {taxid}")
            return None
            
        # Get assembly details
        assembly_id = record["IdList"][0]
        handle = Entrez.esummary(db="assembly", id=assembly_id)
        summary = Entrez.read(handle)
        handle.close()
        
        # Get FTP path - try both RefSeq and GenBank paths
        ftp_path = summary["DocumentSummarySet"]["DocumentSummary"][0].get("FtpPath_RefSeq")
        if not ftp_path:
            ftp_path = summary["DocumentSummarySet"]["DocumentSummary"][0].get("FtpPath_GenBank")
            
        if not ftp_path:
            logger.warning(f"No FTP path found for taxid {taxid}")
            return None
            
        # Construct download URL
        asm_name = os.path.basename(ftp_path)
        url = f"{ftp_path}/{asm_name}_genomic.fna.gz"
        
        # Convert FTP URL to HTTP URL
        http_url = url.replace('ftp://', 'https://')
        
        # Download and save
        try:
            urllib.request.urlretrieve(http_url, genome_path)
            logger.info(f"Successfully downloaded genome for taxid {taxid}")
            return str(genome_path)
        except Exception as e:
            logger.error(f"Failed to download genome for taxid {taxid}: {str(e)}")
            return None
            
    except Exception as e:
        logger.error(f"Error fetching genome for taxid {taxid}: {str(e)}")
        return None