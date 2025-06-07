import unittest
import tempfile
import os
from samovar.fasta_processor import read_fasta, apply_mutations, preprocess_fasta, process_fasta_directories
import pytest
from pathlib import Path
import gzip

class TestFastaProcessor(unittest.TestCase):
    def setUp(self):
        # Create a temporary test FASTA file
        self.test_fasta_content = """>test_seq1
ATCGATCG
>test_seq2
GCTAGCTA"""
        
        self.temp_dir = tempfile.TemporaryDirectory()
        self.test_fasta_path = os.path.join(self.temp_dir.name, 'test.fasta')
        with open(self.test_fasta_path, 'w') as f:
            f.write(self.test_fasta_content)
    
    def tearDown(self):
        self.temp_dir.cleanup()
    
    def test_read_fasta(self):
        sequences = read_fasta(self.test_fasta_path)
        self.assertEqual(len(sequences), 2)
        self.assertEqual(sequences[0][0], 'test_seq1')
        self.assertEqual(sequences[0][1], 'ATCGATCG')
        self.assertEqual(sequences[1][0], 'test_seq2')
        self.assertEqual(sequences[1][1], 'GCTAGCTA')
    
    def test_apply_mutations(self):
        sequence = 'ATCGATCG'
        # Test with mutation rate 0 (should return same sequence)
        self.assertEqual(apply_mutations(sequence, 0.0), sequence)
        
        # Test with mutation rate 1 (should return different sequence)
        mutated = apply_mutations(sequence, 1.0)
        self.assertNotEqual(mutated, sequence)
        self.assertEqual(len(mutated), len(sequence))
        self.assertTrue(all(n in 'ATCG' for n in mutated))
    
    def test_preprocess_fasta(self):
        output_path = os.path.join(self.temp_dir.name, 'output.fasta')
        
        # Test with 50% include and 0 mutation rate
        preprocess_fasta(self.test_fasta_path, output_path, 0.0, 50.0)
        
        with open(output_path, 'r') as f:
            content = f.read()
            
        # Check if output contains expected number of sequences
        # Each input sequence should be split into 2 parts (50% each)
        self.assertEqual(content.count('>test'), 4)
        
        # Check if sequences are properly split
        sequences = content.split('>')[1:]  # Skip empty first element
        for seq in sequences:
            header, sequence = seq.strip().split('\n')
            self.assertEqual(len(sequence), 4)  # Each part should be 4 nucleotides
    
    def test_preprocess_fasta_with_mutations(self):
        output_path = os.path.join(self.temp_dir.name, 'output.fasta')
        
        # Test with 50% include and 1.0 mutation rate
        preprocess_fasta(self.test_fasta_path, output_path, 1.0, 50.0)
        
        with open(output_path, 'r') as f:
            content = f.read()
        
        # Check if all sequences are mutated
        sequences = content.split('>')[1:]
        for seq in sequences:
            header, sequence = seq.strip().split('\n')
            original_part = 'ATCG' if 'test_seq1' in header else 'GCTA'
            # Check that at least one position is different
            self.assertTrue(any(a != b for a, b in zip(sequence, original_part)))
            self.assertEqual(len(sequence), 4)

@pytest.fixture
def test_dir(tmp_path):
    """Create a temporary directory with test FASTA files."""
    # Create test files
    files = {
        "123.fa": ">seq1\nATCG",
        "456.fna": ">seq2\nGCTA",
        "789.fasta": ">seq3\nTAGC",
        "invalid.txt": "not a fasta file",
        "abc.fa": ">seq4\nCGAT"  # non-numeric taxid
    }
    
    # Create gzipped version of one file
    gzipped_content = ">seq5\nATGC"
    with gzip.open(tmp_path / "321.fa.gz", 'wt') as f:
        f.write(gzipped_content)
    
    # Create regular files
    for filename, content in files.items():
        with open(tmp_path / filename, 'w') as f:
            f.write(content)
    
    return tmp_path

def test_process_fasta_directories(test_dir):
    """Test processing of FASTA files from directories."""
    result = process_fasta_directories([str(test_dir)])
    
    # Check that we got the expected files
    expected_files = {
        str(test_dir / "123.fa"): "123",
        str(test_dir / "456.fna"): "456",
        str(test_dir / "789.fasta"): "789",
        str(test_dir / "321.fa"): "321"  # decompressed file
    }
    print (expected_files)
    assert result == expected_files
    
    # Verify that gzipped file was decompressed
    assert not (test_dir / "321.fa.gz").exists()
    assert (test_dir / "321.fa").exists()

def test_process_fasta_directories_nonexistent():
    """Test handling of non-existent directory."""
    result = process_fasta_directories(["/nonexistent/path"])
    assert result == {}

def test_process_fasta_directories_empty(tmp_path):
    """Test handling of empty directory."""
    result = process_fasta_directories([str(tmp_path)])
    assert result == {}

if __name__ == '__main__':
    unittest.main() 