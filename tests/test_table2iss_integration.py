import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from samovar.table2iss import process_annotation_table
import pytest

def create_test_fasta(output_dir, taxid):
    """Create a test FASTA file for a specific taxid."""
    genome_path = os.path.join(output_dir, f"{taxid}.fa")
    seq = Seq("ATCG" * 100)  # 400 bp sequence
    record = SeqRecord(seq, id=f"test_seq_{taxid}", description="")
    SeqIO.write(record, genome_path, "fasta")
    return genome_path

def test_process_annotation_table():
    """Test processing annotation table and generating reads."""
    # Create output directory
    output_dir = "tests_outs/test_process_annotation_table"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test table
    table_path = os.path.join(output_dir, "test_table.csv")
    df = pd.DataFrame({
        "seqID": ["a1", "a2", "a3"],
        "taxID_k1_0": ["0", "562", "0"],
        "taxID_k2_1": ["0", "562", "562"]
    })
    df.to_csv(table_path, index=False)
    
    # Create genome directory
    genome_dir = os.path.join(output_dir, "genomes")
    os.makedirs(genome_dir, exist_ok=True)
    
    # Create test genomes for taxids
    for taxid in ["0", "562"]:
        create_test_fasta(genome_dir, taxid)
    
    # Create output directory for reads
    reads_dir = os.path.join(output_dir, "reads")
    os.makedirs(reads_dir, exist_ok=True)
    
    # Process table
    process_annotation_table(
        table_path=table_path,
        genome_dir=genome_dir,
        output_dir=reads_dir,
        total_amount=1000,
        read_length=150,
        model="hiseq",
        email="test@samovar.com"
    )
    
    # Check if output files exist
    for annotator in ["k1", "k2"]:
        assert os.path.exists(os.path.join(reads_dir, f"test_table_{annotator}_R1.fastq"))
        assert os.path.exists(os.path.join(reads_dir, f"test_table_{annotator}_R2.fastq"))

def test_process_abundance_table_integration():
    """Integration test for process_abundance_table with process_annotation_table."""
    # Create output directory
    output_dir = "tests_outs/test_process_abundance_table_integration"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test annotation table
    table_path = os.path.join(output_dir, "test_annotation.csv")
    df = pd.DataFrame({
        "seqID": ["seq1", "seq2", "seq3"],
        "taxID_k1_0": ["562", "562", "9606"],
        "taxID_k2_1": ["562", "9606", "9606"]
    })
    df.to_csv(table_path, index=False)
    
    # Create genome directory
    genome_dir = os.path.join(output_dir, "genomes")
    os.makedirs(genome_dir, exist_ok=True)
    
    # Create test genomes for taxids
    for taxid in ["562", "9606"]:
        create_test_fasta(genome_dir, taxid)
    
    # Create output directory for reads
    reads_dir = os.path.join(output_dir, "reads")
    os.makedirs(reads_dir, exist_ok=True)
    
    # Process annotation table
    process_annotation_table(
        table_path=table_path,
        genome_dir=genome_dir,
        output_dir=reads_dir,
        total_amount=1000,
        read_length=150,
        model="hiseq",
        email="test@samovar.com"
    )
    
    # Check if output files exist
    for annotator in ["k1", "k2"]:
        assert os.path.exists(os.path.join(reads_dir, f"merged_{annotator}_R1.fastq"))
        assert os.path.exists(os.path.join(reads_dir, f"merged_{annotator}_R2.fastq")) 