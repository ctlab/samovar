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
import subprocess
from .genome_fetcher import fetch_genome

def parse_annotation_table(table_path: str) -> pd.DataFrame:
    """
    Parse annotation table and return DataFrame with taxid and counts for each annotation method.
    
    Args:
        table_path: Path to taxonomy table
        
    Returns:
        DataFrame with columns: taxid, N_k1, N_k2, ... where each N_k column contains
        the count of occurrences for that annotation method
    """
    # Read the input table
    df = pd.read_csv(table_path, sep=",")
    
    # Get all taxid columns (they contain 'taxid' in their name)
    taxid_cols = [col for col in df.columns if 'taxid' in col.lower()]
    
    # Create a new DataFrame to store results
    result = pd.DataFrame()
    
    # For each taxid column, calculate counts
    for col in taxid_cols:
        # Get the annotation name from the column name (e.g., 'k1' from 'taxID_k1')
        ann_name = col.split('_')[-2]
        
        # Count occurrences of each taxid
        taxid_counts = df.groupby(col).size().reset_index()
        taxid_counts.columns = ["taxid", f"N_{ann_name}"]
        
        # Merge with result DataFrame
        if result.empty:
            result = taxid_counts
        else:
            result = result.merge(taxid_counts, on='taxid', how='outer')
    
    # Fill NaN values with 0
    result = result.fillna(0)
    
    # Convert taxid to string
    result['taxid'] = result['taxid'].astype(str)
    
    return result

def generate_reads_genome(
    genome_file: str,
    output_file: str,
    amount: int,
    read_length: int = 150,
    model: str = "hiseq"
) -> None:
    """
    Generate simulated reads from a genome file.
    
    Args:
        genome_file: Path to genome FASTA file
        output_file: Path to output FASTQ file
        amount: Number of reads to generate
        read_length: Length of reads to generate
        model: Model to use for simulation
    """

    cmd = f"""
           iss generate \
            --genomes {genome_file} \
            --model {model} \
            --output {output_file} \
            --n_reads {amount}
    """ 

    subprocess.run(cmd, shell=True)

def generate_reads_metagenome(
    genome_files: List[str],
    output_dir: str,
    amount: List[int],
    read_length: int = 150,
    total_amount: int = None,
    sample_name: str = "merged",
    model: str = "hiseq"
) -> None:
    """
    Generate simulated reads from a metagenome.
    
    Args:
        genome_dir: Directory containing genome files
        output_dir: Path to output FASTQ file
        read_length: Length of reads to generate
        amount: Number of reads to generate
        total_amount: Total number of reads to generate
        sample_name: Name of the sample
        model: Model to use for simulation with ISS
    """
    if total_amount is not None:
        amount = [int(N/sum(amount) * total_amount) for N in amount]
    else:
        total_amount = sum(amount)

    for genome_file, N in zip(genome_files, amount):
        generate_reads_genome(genome_file, os.path.join(output_dir, f"{sample_name}"), N, read_length)
    
    # Merge all reads into a single file
    output_file_R1 = os.path.join(output_dir, f"{sample_name}_R1.fastq")
    output_file_R2 = os.path.join(output_dir, f"{sample_name}_R2.fastq")
    cmd = f"""cat {output_dir}/{sample_name}*R1* >> {output_file_R1}
        cat {output_dir}/{sample_name}*R2* >> {output_file_R2}
    """ 
    subprocess.run(cmd, shell=True)

def regenerate_metagenome(
    genome_files: List[str],
    output_dir: str,
    amount: List[int],
    read_length: int = 150,
    total_amount: int = None,
    sample_name: str = None,
    model: str = "hiseq",
    mode: str = "direct"
) -> None:
    """
    Regenerate metagenome reads from a list of genome files.
    """
    if mode == "direct":
        generate_reads_metagenome(
            genome_files = genome_files, 
            output_dir = output_dir, 
            amount = amount,
            read_length = read_length, 
            total_amount = total_amount, 
            model = model, 
            sample_name = sample_name)
    elif mode == "samovar":
        raise NotImplementedError("Samovar mode not implemented yet")

def process_annotation_table(
    table_path: str,
    genome_dir: str,
    output_dir: str,
    total_amount: int = None,
    email: str = "test@samovar.com",
    reference_only: bool = True,
    model: str = "hiseq",
    read_length: int = 150,
    sample_name: str = None,
    mode: str = "direct"
) -> None:
    """
    Process taxonomy table and generate simulated reads for each taxid.
    
    Args:
        table_path: Path to taxonomy table file
        genome_dir: Directory containing genome files
        output_dir: Directory to write output files
        read_length: Length of reads to generate
        email: Email for NCBI Entrez
        reference_only: If True, only fetch reference genome
        model: Model to use for simulation with ISS
        sample_name: Name of the sample
        total_amount: Total number of reads to generate
        mode: Mode to use for simulation. "Direct" to generate exactly same reads as in the table, "Samovar" to generate reads using Samova.R
    """
    # Read taxonomy table
    annotation_table = parse_annotation_table(table_path)
    if sample_name is None:
        sample_name = os.path.basename(table_path).split(".")[0]

    for taxid in annotation_table['taxid']:
        genome_file = os.path.join(genome_dir, f"{taxid}.fa")
        if not os.path.exists(genome_file):
            fetch_genome(taxid, genome_dir, email, reference_only=True)

    N_cols = [col for col in annotation_table.columns if 'n' in col.lower()]
    for N_annotator in N_cols:
        annotator_name = N_annotator.split('_')[-1]
        amount = annotation_table[N_annotator].tolist()
        taxid = annotation_table['taxid'].tolist()
        genome_files = [os.path.join(genome_dir, f"{taxid}.fa") for taxid in taxid]
        
        # Create output directory for this annotator
        annotator_output_dir = os.path.join(output_dir, annotator_name)
        os.makedirs(annotator_output_dir, exist_ok=True)
        
        regenerate_metagenome(
            genome_files=genome_files,
            output_dir=annotator_output_dir,
            amount=amount,
            read_length=read_length,
            total_amount=total_amount,
            sample_name=sample_name,
            model=model,
            mode=mode
        ) 