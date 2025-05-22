"""
Module for converting tables to simulated reads using ISS-like functionality.
"""

import os
import random
from typing import Dict, List, Tuple
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def read_taxonomy_table(table_path: str) -> Dict[str, int]:
    """
    Read taxonomy table and return dictionary of taxid to amount.
    
    Args:
        table_path: Path to taxonomy table
        
    Returns:
        Dictionary mapping taxid to amount
    """
    df = pd.read_csv(table_path, sep="\t")
    # Convert taxid to string and ensure it's the index
    df["taxid"] = df["taxid"].astype(str)
    df = df.set_index("taxid")
    return df["amount"].to_dict()

def generate_reads(
    genome_file: str,
    output_file: str,
    read_length: int = 150,
    coverage: int = 30,
    amount: int = None
) -> None:
    """
    Generate simulated reads from a genome file.
    
    Args:
        genome_file: Path to genome FASTA file
        output_file: Path to output FASTQ file
        read_length: Length of reads to generate
        coverage: Coverage depth to generate
        amount: Number of reads to generate (overrides coverage)
    """
    # Read genome
    genome = next(SeqIO.parse(genome_file, "fasta"))
    genome_length = len(genome.seq)
    
    # Calculate number of reads to generate
    if amount is None:
        amount = int((genome_length * coverage) / read_length)
    
    # Generate reads
    reads = []
    for i in range(amount):
        # Random start position
        start = random.randint(0, genome_length - read_length)
        read_seq = genome.seq[start:start + read_length]
        
        # Generate random quality scores
        quality = ''.join(chr(random.randint(33, 73)) for _ in range(read_length))
        
        # Create SeqRecord
        read = SeqRecord(
            read_seq,
            id=f"read_{i}",
            description="",
            letter_annotations={"phred_quality": [ord(q) - 33 for q in quality]}
        )
        reads.append(read)
    
    # Write reads to FASTQ file
    SeqIO.write(reads, output_file, "fastq")

def process_table(
    table_path: str,
    genome_dir: str,
    output_dir: str,
    read_length: int = 150,
    coverage: int = 30
) -> None:
    """
    Process taxonomy table and generate simulated reads for each taxid.
    
    Args:
        table_path: Path to taxonomy table file
        genome_dir: Directory containing genome files
        output_dir: Directory to write output files
        read_length: Length of reads to generate
        coverage: Coverage depth to generate
    """
    # Read taxonomy table
    taxid_amounts = read_taxonomy_table(table_path)
    
    # Process each taxid
    for taxid, amount in taxid_amounts.items():
        genome_file = os.path.join(genome_dir, f"{taxid}.fa")
        if not os.path.exists(genome_file):
            print(f"Warning: Genome file not found for taxid {taxid}")
            continue
            
        output_file = os.path.join(output_dir, f"simulated_reads_{taxid}.fastq")
        generate_reads(
            genome_file,
            output_file,
            read_length=read_length,
            coverage=coverage,
            amount=amount
        ) 