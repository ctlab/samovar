"""
Module for converting tables to simulated reads using ISS-like functionality.
"""

import os
import re
import glob
import pandas as pd
import subprocess
from .genome_fetcher import fetch_genome
from typing import List
import yaml
import json

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

def get_genome_file(genome_dir: str, taxid: str) -> str:
    """
    Get the path to a genome file, checking multiple possible extensions.
    
    Args:
        genome_dir: Directory containing genome files
        taxid: Taxonomy ID of the genome
        
    Returns:
        Path to the genome file if found, None otherwise
    """
    extensions = ['.fa', '.fna', '.fasta']
    for ext in extensions:
        genome_file = os.path.join(genome_dir, f"{taxid}{ext}")
        if os.path.exists(genome_file):
            return genome_file
    return None

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
    
    if genome_file is not None and os.path.exists(genome_file):

        cmd = f"""
            iss generate \
                --genomes {genome_file} \
                --model {model} \
                --output {output_file} \
                --n_reads {amount} \
                --fragment-length {read_length}
        """ 

        subprocess.run(cmd, shell=True)

def generate_reads_metagenome(
    genome_files: List[str],
    output_dir: str,
    amount: List[int],
    read_length: int = 150,
    total_amount: int = None,
    sample_name: str = "merged",
    annotator_name: str = None,
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
        annotator_name: Name of the annotator
        sample_name: Name of the sample
        model: Model to use for simulation with ISS
    """
    
    if total_amount is not None:
        amount = [int(N/sum(amount) * total_amount) for N in amount]
    else:
        total_amount = sum(amount)

    if annotator_name is None:
        annotator_name = "any"

    for genome_file, N in zip(genome_files, amount):
        if N > 0:
            generate_reads_genome(
                genome_file, 
                os.path.join(
                    output_dir, 
                    f"{sample_name}_{annotator_name}_{os.path.basename(genome_file).split('.')[0]}"), 
                int(N), read_length, model)
    
    # Merge all reads into a single file
    
    output_file_R1 = os.path.join(output_dir, f"{sample_name}_{annotator_name}_R1.fastq")
    output_file_R2 = os.path.join(output_dir, f"{sample_name}_{annotator_name}_R2.fastq")
    cmd = f"""cat $(ls {output_dir}/{sample_name}_{annotator_name}_*_R1.fastq | grep -v ".iss.tmp.") >> {output_file_R1}
        cat $(ls {output_dir}/{sample_name}_{annotator_name}_*_R2.fastq | grep -v ".iss.tmp.") >> {output_file_R2}
    """ 
    subprocess.run(cmd, shell=True)

    # Cleanup
    for file in glob.glob(os.path.join(output_dir, f"{annotator_name}*")):
        if not os.path.basename(file).startswith("full"):
            os.remove(file)

def regenerate_metagenome(
    genome_files: List[str],
    output_dir: str,
    amount: List[int],
    read_length: int = 150,
    total_amount: int = None,
    sample_name: str = None,
    annotator_name: str = None,
    model: str = "hiseq"
) -> None:
    """
    Regenerate metagenome reads from a list of genome files.
    """
    generate_reads_metagenome(
        genome_files=genome_files,
        output_dir=output_dir,
        amount=amount,
        read_length=read_length,
        total_amount=total_amount,
        model=model,
        sample_name=sample_name,
        annotator_name=annotator_name,
    )

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
        genome_file = get_genome_file(genome_dir, taxid)
        if genome_file is None:
            fetch_genome(taxid, genome_dir, email, reference_only=True)

    N_cols = [col for col in annotation_table.columns if 'n' in col.lower()]
    for N_annotator in N_cols:
        annotator_name = re.search(r'N_(.*?)(?:_[0-9]*)?$', N_annotator).group(1)
        amount = annotation_table[N_annotator].tolist()
        taxid = annotation_table['taxid'].tolist()
        genome_files = [get_genome_file(genome_dir, taxid) for taxid in taxid]
        genome_files = [f for f in genome_files if f is not None]  # Filter out None values
        
        os.makedirs(output_dir, exist_ok=True)
        
        regenerate_metagenome(
            genome_files=genome_files,
            output_dir=output_dir,
            amount=amount,
            read_length=read_length,
            total_amount=total_amount,
            sample_name=sample_name,
            model=model,
            annotator_name=annotator_name
        ) 

def samovar_annotation_regenerate(
    config_samovar: str,
    annotation_dir: str,
    output_dir: str = None
) -> None:
    """
    Regenerate taxonomy tables to a SAMOVAR table.

    Args:
        config_samovar: Path to SAMOVAR config file
        annotation_dir: Path to annotation directory
        output_dir: Path to output directory
    """
    config_samovar_dict = yaml.load(open(config_samovar))
    if output_dir is None:
        output_dir = config_samovar_dict['output_dir']

    try:
        here = os.path.dirname(os.path.abspath(__file__))
        R_path = json.load(open(os.path.join(here, '../build/config.json')))['R_path']
    except:
        R_path = "R"

    subprocess.run(f"{R_path} workflow/annotation_regenerate.R {config_samovar} {annotation_dir} {output_dir}", shell=True)