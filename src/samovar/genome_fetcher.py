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
    genome_path = output_path / f"{taxid}.fa"
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

def parse_taxonomy_table(taxonomy_file: str, output_folder: str, email: str, reference_only: bool = True) -> List[str]:
    """
    Parse taxonomy table and download genomes for all unique taxids.
    
    Args:
        taxonomy_file (str): Path to taxonomy table file
        output_folder (str): Path to output folder for genomes
        email (str): Email for NCBI Entrez
        reference_only (bool): If True, only fetch reference genomes
        
    Returns:
        List[str]: List of paths to downloaded genomes
    """
    try:
        # Read taxonomy table
        df = pd.read_csv(taxonomy_file, sep='\t')
        
        # Find columns with 'taxid' in name
        taxid_columns = [col for col in df.columns if 'taxid' in col.lower()]
        
        if not taxid_columns:
            logger.warning("No columns with 'taxid' found in the taxonomy table")
            return []
            
        # Collect unique taxids
        unique_taxids: Set[str] = set()
        for col in taxid_columns:
            unique_taxids.update(df[col].astype(str).unique())
            
        # Remove any non-numeric taxids
        unique_taxids = {tid for tid in unique_taxids if tid.isdigit()}
        
        logger.info(f"Found {len(unique_taxids)} unique taxids")
        
        # Download genomes
        downloaded_genomes = []
        for taxid in unique_taxids:
            genome_path = fetch_genome(taxid, output_folder, email, reference_only)
            if genome_path:
                downloaded_genomes.append(genome_path)
                
        return downloaded_genomes
        
    except Exception as e:
        logger.error(f"Error parsing taxonomy table: {str(e)}")
        return [] 