#!/usr/bin/env python3

import os
import subprocess
import logging
import tempfile
from pathlib import Path
from typing import List, Dict, Union
from Bio import SeqIO
from Bio.Seq import Seq
import yaml
from .fasta_processor import process_fasta_directories
from .fasta_processor import preprocess_fasta

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_command(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and return the result.
    
    Args:
        cmd (List[str]): List of command and arguments to execute.
        check (bool, optional): If True, raises CalledProcessError if command fails. Defaults to True.
    
    Returns:
        subprocess.CompletedProcess: Object containing the result of the command execution.
    
    Raises:
        subprocess.CalledProcessError: If check=True and command returns non-zero exit code.
    """
    logger.info(f"Running command: {' '.join(cmd)}")
    return subprocess.run(cmd, check=check, text=True)

def get_taxonomy_db(
    db_path: str="kraken_db",
    taxonomy_path: str=None
) -> None:
    """Download or copy the Kraken2 taxonomy database.
    
    This function either downloads the NCBI taxonomy database or copies it from a specified path.
    The taxonomy database is required for Kraken2 to function properly.
    
    Args:
        db_path (str, optional): Path where the database will be stored. Defaults to "kraken_db".
        taxonomy_path (str, optional): Path to existing taxonomy database. If None, downloads from NCBI. Defaults to None.
    
    Note:
        If taxonomy_path is None, the function will download the latest taxonomy database from NCBI.
        The downloaded file will be automatically extracted to the specified db_path.
    """
    os.makedirs(db_path, exist_ok=True)

    if taxonomy_path is None:
        download_cmd =[
            "wget",
            "ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
            "-P",
            db_path
        ]
        run_command(download_cmd)

        extract_cmd = [
            "tar",
            "-xzf",
            db_path + "/taxdump.tar.gz",
            "-C",
            db_path
        ]
        print(extract_cmd)
        run_command(extract_cmd)

        logger.info(f"Taxonomy database downloaded and extracted to {db_path}/taxonomy")
    
    else:
        copy_cmd = [
            "cp",
            "-r",
            taxonomy_path,
            db_path
        ]
        run_command(copy_cmd)

        logger.info(f"Taxonomy database copied to {db_path}/taxonomy")

def process_fasta_kraken2(input_file: str, 
                  taxid: str) -> str:
    """Process FASTA file to add proper taxonomy IDs to headers for Kraken2 compatibility.
    
    This function reads a FASTA file and modifies sequence headers to include taxonomy IDs
    in the format required by Kraken2: "seq{index}|kraken:taxid|{taxid}".
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to be added to sequence headers.
    
    Returns:
        str: Path to the temporary processed FASTA file.
    
    Note:
        The function creates a temporary file that should be cleaned up by the caller.
    """
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    
    for i, record in enumerate(SeqIO.parse(input_file, "fasta")):
        record.id = f"seq{i}|kraken:taxid|{taxid}"

    SeqIO.write(record, temp_fasta, "fasta")
    temp_fasta.close()
    return temp_fasta.name

def add_database_kraken2(
    input_file: str,
    taxid: str,
    db_path: str="kraken_db" 
) -> None:
    """Add sequences to the Kraken2 database library.
    
    This function processes a FASTA file and adds its sequences to the Kraken2 database library
    with proper taxonomy IDs. The sequences must be in FASTA format.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to associate with the sequences.
        db_path (str, optional): Path to the Kraken2 database. Defaults to "kraken_db".
    
    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    # Create database directory if it doesn't exist
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    else:
        temp_fasta = process_fasta_kraken2(input_file, taxid)

    add_cmd = [
        "kraken2-build", 
        "--add-to-library", 
        temp_fasta, 
        "--db", 
        db_path
    ]
    run_command(add_cmd)

    logger.info(f"Database successfully added to {db_path}")


def build_database_kraken2(
    db_path: str="kraken_db", 
    threads: int = 1,
    kmer_len: int = 35,
    minimizer_len: int = 31,
    minimizer_spaces: int = 7,
    skip_maps: bool = True
) -> None:
    """Build a custom Kraken2 database from the library of sequences.
    
    This function builds a Kraken2 database using the sequences that have been added to the library.
    It creates the k-mer database and performs necessary indexing for efficient sequence classification.
    
    Args:
        db_path (str, optional): Path to the database directory. Defaults to "kraken_db".
        threads (int, optional): Number of threads to use for building. Defaults to 1.
        kmer_len (int, optional): Length of k-mers to use. Defaults to 35.
        minimizer_len (int, optional): Length of minimizers to use. Defaults to 31.
        minimizer_spaces (int, optional): Number of spaces in minimizer. Defaults to 7.
        skip_maps (bool, optional): Whether to skip downloading accession maps. Defaults to True.
    
    Note:
        The function will clean up intermediate files after building the database.
        Building a large database can be memory-intensive and time-consuming.
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

def process_fasta_kaiju(input_file: str, 
                       taxid: str,
                       db_path: str = "kaiju_db",
                       protein: bool = False) -> str:
    """Process FASTA file for Kaiju database building.
    
    This function reads a FASTA file and either keeps it as is (if protein=True)
    or translates nucleotide sequences to protein sequences in all three frames (if protein=False).
    It also adds taxonomy IDs to the sequence headers.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to be added to sequence headers.
        db_path (str, optional): Path to the database directory. Defaults to "kaiju_db".
        protein (bool, optional): Whether the input is already protein sequences. Defaults to False.
    
    Returns:
        str: Path to the temporary processed FASTA file.
    
    Note:
        The function creates a temporary file that should be cleaned up by the caller.
        For nucleotide sequences, it generates three frames per sequence.
    """
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    
    if protein:
        # If input is already protein, just copy the file with modified headers
        for record in SeqIO.parse(input_file, "fasta"):
            record.id = f"{record.id}|taxid|{taxid}"
            SeqIO.write(record, temp_fasta, "fasta")
    else:
        # Translate nucleotide sequences to protein in all three frames
        for record in SeqIO.parse(input_file, "fasta"):
            # Trim sequence to multiple of 3 to avoid partial codons
            seq_len = len(record.seq)
            trim_len = (seq_len // 3) * 3
            trimmed_seq = record.seq[:trim_len]
            
            # Create three frames
            for frame in range(3):
                # Get the frame sequence
                frame_seq = trimmed_seq[frame:]
                # Trim to multiple of 3
                frame_trim_len = (len(frame_seq) // 3) * 3
                frame_seq = frame_seq[:frame_trim_len]
                
                # Translate the sequence
                protein_seq = str(frame_seq.translate())
                
                # Create new record with translated sequence and taxid
                new_record = record
                new_record.seq = Seq(protein_seq)
                new_record.id = f"{record.id}_frame{frame+1}|taxid|{taxid}"
                SeqIO.write(new_record, temp_fasta, "fasta")
    
    temp_fasta.close()
    return temp_fasta.name

def add_database_kaiju(
    input_file: str,
    taxid: str,
    db_path: str = "kaiju_db",
    protein: bool = False
) -> None:
    """Add sequences to the Kaiju database library.
    
    This function processes a FASTA file and adds its sequences to the Kaiju database library
    with proper taxonomy IDs. The sequences will be translated to protein if they are not already.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to associate with the sequences.
        db_path (str, optional): Path to the Kaiju database. Defaults to "kaiju_db".
        protein (bool, optional): Whether the input is already protein sequences. Defaults to False.
    
    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    # Create database directory if it doesn't exist
    os.makedirs(db_path, exist_ok=True)
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Process the input file
    processed_file = process_fasta_kaiju(input_file, taxid, db_path, protein)
    
    # Append to the library file
    library_file = os.path.join(db_path, "library.faa")
    with open(processed_file, 'r') as infile, open(library_file, 'a') as outfile:
        outfile.write(infile.read())
    
    # Clean up temporary file
    os.remove(processed_file)
    
    logger.info(f"Sequences successfully added to Kaiju database at {db_path}")

def build_database_kaiju(
    db_path: str = "kaiju_db",
    threads: int = 1,
    protein: bool = False
) -> None:
    """Build a custom Kaiju database from the processed sequences.
    
    This function builds a Kaiju database using the sequences that have been processed
    and added to the database directory. It builds the BWT first, then the FM-index
    for efficient sequence classification.
    
    Args:
        db_path (str, optional): Path to the database directory. Defaults to "kaiju_db".
        threads (int, optional): Number of threads to use for building. Defaults to 1.
        protein (bool, optional): Whether the input is already protein sequences. Defaults to False.
    
    Note:
        Building a large database can be memory-intensive and time-consuming.
    """
    # Create database directory if it doesn't exist
    os.makedirs(db_path, exist_ok=True)
    
    library_file = os.path.join(db_path, "library.faa")
    base_name = os.path.join(db_path, "kaiju_db")
    
    # First build the BWT
    mkbwt_cmd = [
        "kaiju-mkbwt",
        "-o", base_name,
        library_file,
    ]
    run_command(mkbwt_cmd)
    
    # Then build the FM-index
    mkfmi_cmd = [
        "kaiju-mkfmi",
        base_name
    ]
    run_command(mkfmi_cmd)
    
    # Clean up temporary files
    sa_file = f"{base_name}.sa"
    if os.path.exists(sa_file):
        os.remove(sa_file)
    
    bwt_file = f"{base_name}.bwt"
    if os.path.exists(bwt_file):
        os.remove(bwt_file)
    
    logger.info(f"Kaiju database successfully built at {db_path}")

def process_fasta_krakenunique(input_file: str, taxid: str, genome_name: str, db_path: str) -> tuple[str, str]:
    """Process FASTA file for KrakenUniq database building.
    
    This function preprocesses a FASTA file and creates a genomes.map file entry.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to be added to sequence headers.
        genome_name (str): Name of the genome/assembly.
        db_path (str): Path to the database directory.
    
    Returns:
        tuple[str, str]: Paths to the processed FASTA file and genomes.map file.
    """
    # Preprocess the FASTA file
    processed_fasta = os.path.join(db_path, os.path.basename(input_file))
    preprocess_fasta(input_file, processed_fasta, mutation_rate=0, include_percent=100)
    
    # Create genomes.map entry
    map_file = os.path.join(db_path, "genomes.map")
    with open(map_file, 'a') as map_out:
        for record in SeqIO.parse(processed_fasta, "fasta"):
            # Get sequence ID (header without '>' up to first space)
            seq_id = record.id.split()[0]
            # Write to genomes.map
            map_out.write(f"{seq_id}\t{taxid}\t{genome_name}\n")
    
    return processed_fasta, map_file

def add_database_krakenunique(
    input_file: str,
    taxid: str,
    genome_name: str,
    db_path: str = "krakenuniq_db"
) -> None:
    """Add sequences to the KrakenUniq database library.
    
    This function processes a FASTA file and adds its sequences to the KrakenUniq database library
    with proper taxonomy IDs and creates a genomes.map file.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to associate with the sequences.
        genome_name (str): Name of the genome/assembly.
        db_path (str, optional): Path to the KrakenUniq database. Defaults to "krakenuniq_db".
    
    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    # Create database directory if it doesn't exist
    os.makedirs(db_path, exist_ok=True)
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Process the input file
    process_fasta_krakenunique(input_file, taxid, genome_name, db_path)
    
    logger.info(f"Sequences successfully added to KrakenUniq database at {db_path}")

def build_database_krakenunique(
    db_path: str = "krakenuniq_db",
    threads: int = 1
) -> None:
    """Build a custom KrakenUniq database from the library of sequences.
    
    This function builds a KrakenUniq database using the sequences that have been added to the library.
    
    Args:
        db_path (str, optional): Path to the database directory. Defaults to "krakenuniq_db".
        threads (int, optional): Number of threads to use for building. Defaults to 1.
    
    Note:
        Building a large database can be memory-intensive and time-consuming.
    """
    # Check if krakenuniq is installed
    try:
        subprocess.run(['which', 'krakenuniq-build'], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        raise RuntimeError("krakenuniq-build not found. Please install KrakenUniq first.")

    # Check if jellyfish is installed
    try:
        jellyfish_path = subprocess.run(['which', 'jellyfish'], check=True, capture_output=True, text=True).stdout.strip()
    except subprocess.CalledProcessError:
        raise RuntimeError("jellyfish not found. Please install Jellyfish first.")

    # Set up environment variables
    env = os.environ.copy()
    env['JELLYFISH_BIN'] = jellyfish_path

    build_cmd = [
        "krakenuniq-build",
        "--db", db_path,
        "--threads", str(threads),
        "--taxids-for-genomes",
        "--taxids-for-sequences"
    ]
    
    # Run command with modified environment
    logger.info(f"Running command: {' '.join(build_cmd)}")
    result = subprocess.run(build_cmd, check=True, text=True, env=env)
    
    logger.info(f"KrakenUniq database successfully built at {db_path}")

def build_database_from_config(config_path: str, db_type: str = "kaiju", db_path: str = "tests_outs/db"):
    """
    Build database from config file.
    
    Args:
        config_path: Path to config YAML file
        db_type: Type of database to build ("kaiju", "kraken2", or "krakenunique")
        db_path: Path to store the database
    """
    # Load configuration
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Process input directories
    file_taxid_map = process_fasta_directories(config['input_dir'])
    
    # Get lists of files and taxids
    input_files = list(file_taxid_map.keys())
    taxids = list(file_taxid_map.values())
    
    print(f"Processing {len(input_files)} files with taxids: {taxids}")
    
    # Build database based on type
    if db_type == "kaiju":
        for input_file, taxid in zip(input_files, taxids):
            add_database_kaiju(input_file, taxid, db_path=db_path)
        get_taxonomy_db(db_path=db_path)
        build_database_kaiju(db_path=db_path, threads=1, protein=False)
    elif db_type == "kraken2":
        for input_file, taxid in zip(input_files, taxids):
            add_database_kraken2(input_file, taxid, db_path=db_path)
        get_taxonomy_db(db_path=db_path+"/taxonomy")
        build_database_kraken2(db_path=db_path, threads=1, kmer_len=35, 
                             minimizer_len=31, minimizer_spaces=7, skip_maps=True)
    elif db_type == "krakenunique":
        # Create output directory
        os.makedirs(db_path, exist_ok=True)
        # Get taxonomy database
        get_taxonomy_db(db_path=db_path)
        # Add each genome to the database
        for input_file, taxid in zip(input_files, taxids):
            # Use filename as genome name
            genome_name = os.path.basename(input_file)
            add_database_krakenunique(input_file, taxid, genome_name, db_path=db_path)
        # Build the database
        build_database_krakenunique(db_path=db_path, threads=1)
    else:
        raise ValueError(f"Unknown database type: {db_type}")