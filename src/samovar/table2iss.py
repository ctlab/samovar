"""
Module for converting tables to simulated reads using ISS-like functionality.
"""

import os
import re
import glob
import pandas as pd
import subprocess
import warnings
import yaml
import json
import tempfile
from typing import List, Union
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
    df = pd.read_csv(table_path, sep="\t")
    
    # Get all taxid columns (they contain 'taxid' in their name)
    taxid_cols = [col for col in df.columns if 'taxid' in col.lower()]
    
    # Create a new DataFrame to store results
    result = pd.DataFrame()
    
    # For each taxid column, calculate counts
    for col in taxid_cols:
        ann_name = col.lower().replace('taxid_', '')
        
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
    extensions = ['-processed.fa', '-processed.fna', '-processed.fasta',
                  '.fa', '.fna', '.fasta', '.fa.gz', '.fna.gz', '.fasta.gz']
    
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
    """
    if genome_file is not None and os.path.exists(genome_file):
        cmd = f"""
            iss generate \
                --genomes {genome_file} \
                --model {model} \
                --output {output_file} \
                --n_reads {amount} \
                --cpus 1
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
    """
    if total_amount is not None:
        amount = [int(N/sum(amount) * total_amount) for N in amount]
    else:
        total_amount = sum(amount)

    if annotator_name is None:
        annotator_name = "any"

    for genome_file, N in zip(genome_files, amount):
        if N > 0:
            output_prefix = os.path.join(
                output_dir, 
                f"{sample_name}_{annotator_name}_{os.path.basename(genome_file).split('.')[0]}"
            )
            
            if N <= 2:
                with open(f"{output_prefix}_R1.fastq", "w") as f1, open(f"{output_prefix}_R2.fastq", "w") as f2:
                    f1.write("")
                    f2.write("")
            else:
                generate_reads_genome(genome_file, output_prefix, int(N), read_length, model)
    
    # Merge all reads into a single file
    output_file_R1 = os.path.join(output_dir, f"{sample_name}_{annotator_name}_R1.fastq")
    output_file_R2 = os.path.join(output_dir, f"{sample_name}_{annotator_name}_R2.fastq")
    valid_r1 = [f for f in glob.glob(os.path.join(output_dir, f"{sample_name}_{annotator_name}_*_R1.fastq")) if ".iss.tmp." not in f]
    valid_r2 = [f for f in glob.glob(os.path.join(output_dir, f"{sample_name}_{annotator_name}_*_R2.fastq")) if ".iss.tmp." not in f]

    if not valid_r1:
        subprocess.run(f'echo "" > {output_file_R1}', shell=True)
        subprocess.run(f'echo "" > {output_file_R2}', shell=True)
    else:
        subprocess.run(f'cat {" ".join(valid_r1)} >> {output_file_R1}', shell=True)
        subprocess.run(f'cat {" ".join(valid_r2)} >> {output_file_R2}', shell=True)

    # Cleanup (Removed the duplicated buggy block)
    for file in glob.glob(os.path.join(output_dir, f"{annotator_name}*")):
        if not os.path.basename(file).startswith("full"):
            try:
                os.remove(file)
            except OSError:
                pass

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
    sample_name: str = None
) -> None:
    """
    Process taxonomy table and generate simulated reads for each taxid. Inherits from process_abundance_table.
    """
    abundance_table = parse_annotation_table(table_path)
    process_abundance_table(
        table=abundance_table,
        genome_dir=genome_dir,
        output_dir=output_dir,
        total_amount=total_amount,
        email=email,
        reference_only=reference_only,
        model=model,
        read_length=read_length,
        sample_name=sample_name
    )
        
def process_abundance_table(
    table: Union[str, pd.DataFrame],
    genome_dir: str,
    output_dir: str,
    total_amount: int = None,
    email: str = "test@samovar.com",
    reference_only: bool = True,
    model: str = "hiseq",
    read_length: int = 150,
    sample_name: str = None
) -> pd.DataFrame:
    """
    Process abundance table and generate simulated reads for each taxid.
    """
    if isinstance(table, str):
        abundance_table = pd.read_csv(table, sep="\t")
    else:
        abundance_table = table
    
    if sample_name is None:
        if isinstance(table, str):
            sample_name = os.path.basename(table).split(".")[0]
        else:
            sample_name = "merged"

    # Try to fetch missing genomes
    available_genomes = []
    print(f"DEBUG: Processing taxids from table: {abundance_table['taxid'].unique()}")
    for taxid in abundance_table['taxid']:
        taxid = str(taxid).split(".")[0]
        genome_file = get_genome_file(genome_dir, taxid)
        if genome_file is None:
            print(f"DEBUG: Genome {taxid} not found in {genome_dir}, fetching...")
            genome_file = fetch_genome(taxid, genome_dir, email, reference_only=True)
        if genome_file is not None:
            print(f"DEBUG: Found genome for {taxid}: {genome_file}")
            available_genomes.append((taxid, genome_file))
        else:
            print(f"DEBUG: Still could not find/fetch genome for {taxid}")

    # [FIX]: Use warnings.warn instead of raise Warning so it actually returns None
    if not available_genomes:
        print("DEBUG: available_genomes is empty!")
        os.makedirs(output_dir, exist_ok=True)
        file_r1 = os.path.join(output_dir, f"{sample_name}_R1.fastq")
        with open(file_r1, "w") as f:
            f.write("")
        file_r2 = os.path.join(output_dir, f"{sample_name}_R2.fastq")
        with open(file_r2, "w") as f:
            f.write("")
        
        warnings.warn("No genome files available for any taxid")
        return None

    # Filter abundance table to only include available genomes
    available_taxids = [str(taxid).split(".")[0] for taxid, _ in available_genomes]
    abundance_table['taxid'] = abundance_table['taxid'].astype(str).str.split(".").str[0]
    filtered_table = abundance_table[abundance_table['taxid'].isin(available_taxids)]
    
    N_cols = [col for col in filtered_table.columns if 'n' in col.lower()]
    for N_annotator in N_cols:
        annotator_name = re.search(r'N_(.*?)(?:_[0-9]*)?$', N_annotator).group(1)
        amount = filtered_table[N_annotator].tolist()
        taxid = filtered_table['taxid'].tolist()
        genome_files = [get_genome_file(genome_dir, t) for t in taxid]
        genome_files = [f for f in genome_files if f is not None]
        
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
    
    return filtered_table

def samovar_annotation_regenerate(
    annotation_dir: str,
    config_samovar: str = None,
    output_dir: str = None
) -> None:
    """
    Regenerate taxonomy tables to a SAMOVAR table.
    """
    if config_samovar is None:
        tmp_file = tempfile.mktemp()
        with open(tmp_file, 'w') as f:
            yaml.dump({
                'threshold_amount': 1e-5,
                'plot_log': False,
                'N': 10,
                'N_reads': 1000
            }, f)
        config_samovar = tmp_file

    with open(config_samovar, 'r') as f:
        config_samovar_dict = yaml.safe_load(f)

    if output_dir is None:
        output_dir = config_samovar_dict['output_dir']

    try:
        here = os.path.dirname(os.path.abspath(__file__))
        R_path = json.load(open(os.path.join(here, '../../build/config.json')))['r_path']
    except:
        R_path = "R"

    annotation_regenerate = os.path.join(here, '../../workflow/annotation_regenerate.R')
    cmd = f"{R_path} \
        --vanilla -s -f {annotation_regenerate} \
        --args \
        --config {config_samovar} \
        --annotation_dir {annotation_dir} \
        --output_dir {output_dir}"

    subprocess.run(cmd, shell=True)