import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from samovar.table2iss import (
    generate_reads_genome,
    generate_reads_metagenome,
    regenerate_metagenome
)

def create_test_fasta(output_dir):
    """Create a test FASTA file."""
    fasta_path = os.path.join(output_dir, "test.fa")
    seq = Seq("ATCG" * 100)  # 400 bp sequence
    record = SeqRecord(seq, id="test_seq", description="")
    SeqIO.write(record, fasta_path, "fasta")
    return fasta_path

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
        sample_name="test_regenerated"
    )
    
    # Check if output files exist
    assert os.path.exists(os.path.join(output_dir, "test_regenerated_any_R1.fastq"))
    assert os.path.exists(os.path.join(output_dir, "test_regenerated_any_R2.fastq")) 