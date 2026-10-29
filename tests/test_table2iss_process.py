import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from samovar.table2iss import _iss_worker_count, process_abundance_table
import pytest


@pytest.mark.mandatory
def test_iss_workers_match_chunks(tmp_path):
    # 38 reads → 19 pairs. 24 CPUs would leave empty temp files.
    assert _iss_worker_count(24, 38) == 19
    assert _iss_worker_count(12, 38) == 10
    assert _iss_worker_count(4, 1000) == 4
    assert _iss_worker_count(8, 0) == 1
    assert _iss_worker_count(8, 2) == 1

def create_test_fasta(output_dir, taxid):
    """Create a test FASTA file for a specific taxid."""
    genome_path = os.path.join(output_dir, f"{taxid}.fa")
    seq = Seq("ATCG" * 100)  # 400 bp sequence
    record = SeqRecord(seq, id=f"test_seq_{taxid}", description="")
    SeqIO.write(record, genome_path, "fasta")
    return genome_path

@pytest.mark.optional
def test_process_abundance_table_with_dataframe(tmp_path):
    """Test process_abundance_table with DataFrame input."""
    # Create output directory
    output_dir = str(tmp_path / "run")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test DataFrame
    df = pd.DataFrame({
        "taxid": ["562", "562", "9606"],
        "N_k1": [100, 200, 300],
        "N_k2": [150, 250, 350]
    })
    
    # Create genome directory
    genome_dir = os.path.join(output_dir, "genomes")
    os.makedirs(genome_dir, exist_ok=True)
    
    # Create test genomes for taxids
    for taxid in ["562", "9606"]:
        create_test_fasta(genome_dir, taxid)
    
    # Create output directory for reads
    reads_dir = os.path.join(output_dir, "reads")
    os.makedirs(reads_dir, exist_ok=True)
    
    # Process table
    process_abundance_table(
        table=df,
        genome_dir=genome_dir,
        output_dir=reads_dir,
        total_amount=40,
        read_length=150,
        model="hiseq",
        email="test@samovar.com"
    )
    
    # Check if output files exist
    for annotator in ["k1", "k2"]:
        assert os.path.exists(os.path.join(reads_dir, f"merged_{annotator}_R1.fastq"))
        assert os.path.exists(os.path.join(reads_dir, f"merged_{annotator}_R2.fastq"))

@pytest.mark.optional
def test_process_abundance_table_with_file(tmp_path):
    """Test process_abundance_table with file input."""
    # Create output directory
    output_dir = str(tmp_path / "run")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test table
    table_path = os.path.join(output_dir, "test_table.csv")
    df = pd.DataFrame({
        "taxid": ["562", "562", "9606"],
        "N_k1": [100, 200, 300],
        "N_k2": [150, 250, 350]
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
    
    # Process table
    process_abundance_table(
        table=table_path,
        genome_dir=genome_dir,
        output_dir=reads_dir,
        total_amount=40,
        read_length=150,
        model="hiseq",
        email="test@samovar.com",
        sample_name="test_sample"
    )
    
    # Check if output files exist
    for annotator in ["k1", "k2"]:
        assert os.path.exists(os.path.join(reads_dir, f"test_sample_{annotator}_R1.fastq"))
        assert os.path.exists(os.path.join(reads_dir, f"test_sample_{annotator}_R2.fastq"))

@pytest.mark.optional
def test_process_abundance_table_with_missing_genome(tmp_path):
    """Test process_abundance_table with some missing genome files."""
    # Create output directory
    output_dir = str(tmp_path / "run")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test DataFrame
    df = pd.DataFrame({
        "taxid": ["562", "562", "999999999"],
        "N_k1": [100, 200, 300],
        "N_k2": [150, 250, 350]
    })
    
    # Create genome directory
    genome_dir = os.path.join(output_dir, "genomes")
    os.makedirs(genome_dir, exist_ok=True)
    
    # Create test genome for only one taxid
    create_test_fasta(genome_dir, "562")
    
    # Create output directory for reads
    reads_dir = os.path.join(output_dir, "reads")
    os.makedirs(reads_dir, exist_ok=True)
    
    # Process table - should work with at least one genome
    process_abundance_table(
        table=df,
        genome_dir=genome_dir,
        output_dir=reads_dir,
        total_amount=40,
        read_length=150,
        model="hiseq",
        email="test@samovar.com"
    )
    
    # Check if output files exist
    for annotator in ["k1", "k2"]:
        assert os.path.exists(os.path.join(reads_dir, f"merged_{annotator}_R1.fastq"))
        assert os.path.exists(os.path.join(reads_dir, f"merged_{annotator}_R2.fastq"))

@pytest.mark.optional
def test_process_abundance_table_with_zero_reads(tmp_path):
    """Test process_abundance_table with zero reads for some taxids."""
    # Create output directory
    output_dir = str(tmp_path / "run")
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test DataFrame with zero reads for some taxids
    df = pd.DataFrame({
        "taxid": ["562", "562", "9606"],
        "N_k1": [0, 200, 0],
        "N_k2": [150, 0, 350]
    })
    
    # Create genome directory
    genome_dir = os.path.join(output_dir, "genomes")
    os.makedirs(genome_dir, exist_ok=True)
    
    # Create test genomes for taxids
    for taxid in ["562", "9606"]:
        create_test_fasta(genome_dir, taxid)
    
    # Create output directory for reads
    reads_dir = os.path.join(output_dir, "reads")
    os.makedirs(reads_dir, exist_ok=True)
    
    # Process table
    process_abundance_table(
        table=df,
        genome_dir=genome_dir,
        output_dir=reads_dir,
        total_amount=40,
        read_length=150,
        model="hiseq",
        email="test@samovar.com"
    )
    
    # Check if output files exist
    for annotator in ["k1", "k2"]:
        assert os.path.exists(os.path.join(reads_dir, f"merged_{annotator}_R1.fastq"))
        assert os.path.exists(os.path.join(reads_dir, f"merged_{annotator}_R2.fastq")) 