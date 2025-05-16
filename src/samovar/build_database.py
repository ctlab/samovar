#!/usr/bin/env python3

import argparse
import os
import subprocess
import logging
import tempfile
from pathlib import Path
from typing import List, Dict, Union
from Bio import SeqIO

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_command(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and return the result."""
    logger.info(f"Running command: {' '.join(cmd)}")
    return subprocess.run(cmd, check=check, text=True)

def get_taxonomy_db(
    db_path: str="kraken_db",
    taxonomy_path: str=None
) -> None:
    """Get the taxonomy database."""
    os.makedirs(db_path + "/taxonomy", exist_ok=True)

    if taxonomy_path is None:
        download_cmd =[
            "wget",
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            "-P",
            db_path + "/taxonomy"
        ]
        run_command(download_cmd)

        extract_cmd = [
            "tar",
            "-xzf",
            db_path + "/taxonomy/taxdump.tar.gz",
            "-C",
            db_path + "/taxonomy"
        ]
        run_command(extract_cmd)

        logger.info(f"Taxonomy database downloaded and extracted to {db_path}/taxonomy")
    
    else:
        copy_cmd = [
            "cp",
            "-r",
            taxonomy_path,
            db_path + "/taxonomy"
        ]
        run_command(copy_cmd)

        logger.info(f"Taxonomy database copied to {db_path}/taxonomy")

def process_fasta(input_file: str, 
                  taxid: str) -> None:
    """Process FASTA file to add proper taxonomy IDs to headers."""
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    
    for i, record in enumerate(SeqIO.parse(input_file, "fasta")):
        record.id = f"seq{i}|kraken:taxid|{taxid}"

    SeqIO.write(record, temp_fasta, "fasta")
    temp_fasta.close()
    return temp_fasta.name

def add_database(
    input_file: str,
    taxid: str,
    db_path: str="kraken_db" 
) -> None:
    """
    Add input files to the kraken2 database.
    """
    # Create database directory if it doesn't exist
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    else:
        temp_fasta = process_fasta(input_file, taxid)

    add_cmd = [
        "kraken2-build", 
        "--add-to-library", 
        temp_fasta, 
        "--db", 
        db_path
    ]
    run_command(add_cmd)

    logger.info(f"Database successfully added to {db_path}")

    
def build_database(
    db_path: str="kraken_db", 
    threads: int = 1,
    kmer_len: int = 35,
    minimizer_len: int = 31,
    minimizer_spaces: int = 7,
    skip_maps: bool = True
) -> None:
    """
    Build a custom Kraken2 database from input files.
    
    Args:
        db_path: Path to the database directory
        threads: Number of threads to use for building
        kmer_len: Length of k-mers to use
        minimizer_len: Length of minimizers to use
        minimizer_spaces: Number of spaces in minimizer
        skip_maps: Whether to skip downloading accession maps. Not implemented yet.
    """
    
    # Build the database
    build_cmd = [
        "kraken2-build",
        "--build",
        "--db", db_path,
        "--threads", str(threads),
        "--kmer-len", str(kmer_len),
        "--minimizer-len", str(minimizer_len),
        "--minimizer-spaces", str(minimizer_spaces),
        "--skip-maps"
    ]
    run_command(build_cmd)
    
    # Clean up intermediate files
    clean_cmd = ["kraken2-build", "--clean", "--db", db_path]
    run_command(clean_cmd)
    
    logger.info(f"Database successfully built at {db_path}")