import os
import tempfile
import pandas as pd
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from samovar.table2iss import read_taxonomy_table, generate_reads, process_table

def create_test_fasta(tmpdir):
    """Create a test FASTA file."""
    fasta_path = os.path.join(tmpdir, "test.fa")
    seq = Seq("ATCG" * 100)  # 400 bp sequence
    record = SeqRecord(seq, id="test_seq", description="")
    SeqIO.write(record, fasta_path, "fasta")
    return fasta_path

def create_test_table(tmpdir):
    """Create a test taxonomy table."""
    table_path = os.path.join(tmpdir, "test_table.tsv")
    df = pd.DataFrame({
        "taxid": ["9606", "10090", "10116"],
        "amount": [1000, 2000, 3000]
    })
    df.to_csv(table_path, sep="\t", index=False)
    return table_path

def test_read_taxonomy_table():
    """Test reading taxonomy table."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test table
        table_path = create_test_table(tmpdir)
        
        # Read table
        taxid_amounts = read_taxonomy_table(table_path)
        
        # Check results
        assert isinstance(taxid_amounts, dict)
        assert len(taxid_amounts) == 3
        assert taxid_amounts["9606"] == 1000
        assert taxid_amounts["10090"] == 2000
        assert taxid_amounts["10116"] == 3000

def test_generate_reads():
    """Test generating reads from a genome."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test FASTA
        fasta_path = create_test_fasta(tmpdir)
        output_path = os.path.join(tmpdir, "test_reads.fastq")
        
        # Generate reads
        generate_reads(
            fasta_path,
            output_path,
            read_length=100,
            coverage=10
        )
        
        # Check if output file exists
        assert os.path.exists(output_path)
        
        # Read generated reads
        reads = list(SeqIO.parse(output_path, "fastq"))
        
        # Check number of reads (should be approximately genome_length * coverage / read_length)
        expected_reads = (400 * 10) // 100  # 40 reads
        assert len(reads) == expected_reads
        
        # Check read length
        for read in reads:
            assert len(read.seq) == 100
            assert len(read.letter_annotations["phred_quality"]) == 100

def test_generate_reads_with_amount():
    """Test generating reads with specific amount."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test FASTA
        fasta_path = create_test_fasta(tmpdir)
        output_path = os.path.join(tmpdir, "test_reads.fastq")
        
        # Generate reads with specific amount
        generate_reads(
            fasta_path,
            output_path,
            read_length=100,
            amount=50
        )
        
        # Read generated reads
        reads = list(SeqIO.parse(output_path, "fastq"))
        
        # Check number of reads
        assert len(reads) == 50

def test_process_table():
    """Test processing taxonomy table and generating reads."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test table
        table_path = create_test_table(tmpdir)
        
        # Create test genome directory
        genome_dir = os.path.join(tmpdir, "genomes")
        os.makedirs(genome_dir)
        
        # Create test genomes
        for taxid in ["9606", "10090", "10116"]:
            fasta_path = os.path.join(genome_dir, f"{taxid}.fa")
            seq = Seq("ATCG" * 100)
            record = SeqRecord(seq, id=f"seq_{taxid}", description="")
            SeqIO.write(record, fasta_path, "fasta")
        
        # Create output directory
        output_dir = os.path.join(tmpdir, "output")
        os.makedirs(output_dir)
        
        # Process table
        process_table(
            table_path,
            genome_dir,
            output_dir,
            read_length=100,
            coverage=10
        )
        
        # Check output files
        for taxid in ["9606", "10090", "10116"]:
            output_file = os.path.join(output_dir, f"simulated_reads_{taxid}.fastq")
            assert os.path.exists(output_file)
            
            # Check number of reads
            reads = list(SeqIO.parse(output_file, "fastq"))
            # Read table and convert taxid to string for comparison
            table_df = pd.read_csv(table_path, sep="\t")
            table_df["taxid"] = table_df["taxid"].astype(str)
            expected_reads = int(table_df.set_index("taxid").loc[taxid, "amount"])
            assert len(reads) == expected_reads

def test_process_table_missing_genome():
    """Test processing table with missing genome file."""
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test table
        table_path = create_test_table(tmpdir)
        
        # Create test genome directory
        genome_dir = os.path.join(tmpdir, "genomes")
        os.makedirs(genome_dir)
        
        # Create only one test genome
        fasta_path = os.path.join(genome_dir, "9606.fa")
        seq = Seq("ATCG" * 100)
        record = SeqRecord(seq, id="seq_9606", description="")
        SeqIO.write(record, fasta_path, "fasta")
        
        # Create output directory
        output_dir = os.path.join(tmpdir, "output")
        os.makedirs(output_dir)
        
        # Process table
        process_table(
            table_path,
            genome_dir,
            output_dir,
            read_length=100,
            coverage=10
        )
        
        # Check that only one output file exists
        assert os.path.exists(os.path.join(output_dir, "simulated_reads_9606.fastq"))
        assert not os.path.exists(os.path.join(output_dir, "simulated_reads_10090.fastq"))
        assert not os.path.exists(os.path.join(output_dir, "simulated_reads_10116.fastq")) 