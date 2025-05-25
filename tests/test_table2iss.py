import os
import tempfile
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from samovar.table2iss import (
    parse_annotation_table,
    generate_reads_genome,
    generate_reads_metagenome,
    regenerate_metagenome,
    process_annotation_table
)

def create_test_fasta(output_dir):
    """Create a test FASTA file."""
    fasta_path = os.path.join(output_dir, "test.fa")
    seq = Seq("ATCG" * 100)  # 400 bp sequence
    record = SeqRecord(seq, id="test_seq", description="")
    SeqIO.write(record, fasta_path, "fasta")
    return fasta_path

def create_test_table(output_dir):
    """Create a test taxonomy table."""
    table_path = os.path.join(output_dir, "test_table.tsv")
    df = pd.DataFrame({
        "taxid": ["9606", "10090", "10116"],
        "amount": [1000, 2000, 3000]
    })
    df.to_csv(table_path, sep="\t", index=False)
    return table_path

def test_parse_annotation_table():
    """Test reading taxonomy table and counting occurrences."""
    # Create output directory
    output_dir = "tests_outs/test_parse_annotation_table"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test table with the example data
    table_path = os.path.join(output_dir, "test_table.csv")
    df = pd.DataFrame({
        "seqID": ["a1", "a2", "a3"],
        "taxID_k1_0": ["0", "562", "0"],
        "taxID_k2_1": ["0", "562", "562"]
    })
    df.to_csv(table_path, index=False)
    
    # Read table
    result = parse_annotation_table(table_path)
    
    # Save output for inspection
    result.to_csv(os.path.join(output_dir, "result.csv"), index=False)
    
    # Check results
    assert isinstance(result, pd.DataFrame)
    assert set(result.columns) == {'taxid', 'N_k1', 'N_k2'}
    assert len(result) == 2  # 0 and 562
    
    # Check counts
    taxid_0 = result[result['taxid'] == '0'].iloc[0]
    assert taxid_0['N_k1'] == 2
    assert taxid_0['N_k2'] == 1 
    
    taxid_562 = result[result['taxid'] == '562'].iloc[0]
    assert taxid_562['N_k1'] == 1
    assert taxid_562['N_k2'] == 2 

def test_parse_real_annotation_table():
    """Test parsing the real annotation table and check output dimensions."""
    # Create output directory
    output_dir = "tests_outs/test_parse_real_annotation_table"
    os.makedirs(output_dir, exist_ok=True)
    
    # Path to real annotation file
    table_path = "tests/data/annotation.csv"
    
    # Parse table
    result = parse_annotation_table(table_path)
    
    # Save output for manual inspection
    output_path = os.path.join(output_dir, "result.csv")
    result.to_csv(output_path, index=False)
    
    # Check basic properties
    assert isinstance(result, pd.DataFrame)
    assert 'taxid' in result.columns
    assert all(col.startswith('N_') for col in result.columns if col != 'taxid')
    assert len(result) > 0  # Should have at least one taxid

def test_generate_reads_genome():
    """Test generating reads from a single genome."""
    # Create output directory
    output_dir = "tests_outs/test_generate_reads_genome"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create test FASTA file
    fasta_path = create_test_fasta(output_dir)
    
    # Set up output path
    output_path = os.path.join(output_dir, "test_reads")
    
    # Generate reads
    generate_reads_genome(
        genome_file=fasta_path,
        output_file=output_path,
        amount=100,
        read_length=150
    )
    
    # Check if output files exist
    assert os.path.exists(f"{output_path}_R1.fastq")
    assert os.path.exists(f"{output_path}_R2.fastq")

def test_generate_reads_metagenome():
    """Test generating reads from multiple genomes."""
    # Create output directory
    output_dir = "tests_outs/test_generate_reads_metagenome"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create multiple test FASTA files
    fasta_paths = []
    for i in range(3):
        fasta_path = create_test_fasta(output_dir)
        fasta_paths.append(fasta_path)
    
    # Generate reads
    generate_reads_metagenome(
        genome_files=fasta_paths,
        output_dir=output_dir,
        amount=[100, 200, 300],
        read_length=150,
        sample_name="test_metagenome"
    )
    
    # Check if output files exist
    assert os.path.exists(os.path.join(output_dir, "test_metagenome_any_R1.fastq"))
    assert os.path.exists(os.path.join(output_dir, "test_metagenome_any_R2.fastq"))

def test_regenerate_metagenome():
    """Test regenerating metagenome reads."""
    # Create output directory
    output_dir = "tests_outs/test_regenerate_metagenome"
    os.makedirs(output_dir, exist_ok=True)
    
    # Create multiple test FASTA files
    fasta_paths = []
    for i in range(3):
        fasta_path = create_test_fasta(output_dir)
        fasta_paths.append(fasta_path)
    
    # Regenerate reads
    regenerate_metagenome(
        genome_files=fasta_paths,
        output_dir=output_dir,
        amount=[100, 200, 300],
        read_length=150,
        sample_name="test_regenerated",
        mode="direct"
    )
    
    # Check if output files exist
    assert os.path.exists(os.path.join(output_dir, "test_regenerated_any_R1.fastq"))
    assert os.path.exists(os.path.join(output_dir, "test_regenerated_any_R2.fastq"))

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
        genome_path = os.path.join(genome_dir, f"{taxid}.fa")
        seq = Seq("ATCG" * 100)  # 400 bp sequence
        record = SeqRecord(seq, id=f"test_seq_{taxid}", description="")
        SeqIO.write(record, genome_path, "fasta")
    
    # Create output directory for reads
    reads_dir = os.path.join(output_dir, "reads")
    os.makedirs(reads_dir, exist_ok=True)
    
    # Process table
    process_annotation_table(
        table_path=table_path,
        genome_dir=genome_dir,
        output_dir=reads_dir,
        total_amount=1000,
        mode="direct"
    )
    
    # Check if output files exist
    for annotator in ["k1", "k2"]:
        assert os.path.exists(os.path.join(reads_dir, f"test_table_{annotator}_R1.fastq"))
        assert os.path.exists(os.path.join(reads_dir, f"test_table_{annotator}_R2.fastq"))


