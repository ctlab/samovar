import unittest
import tempfile
import os
from pathlib import Path
from samovar.fasta_processor import read_fasta, apply_mutations, process_fasta

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
    
    def test_process_fasta(self):
        output_path = os.path.join(self.temp_dir.name, 'output.fasta')
        
        # Test with 50% include and 0 mutation rate
        process_fasta(self.test_fasta_path, output_path, 0.0, 50.0)
        
        with open(output_path, 'r') as f:
            content = f.read()
            
        # Check if output contains expected number of sequences
        # Each input sequence should be split into 2 parts (50% each)
        self.assertEqual(content.count('>test-'), 4)
        
        # Check if sequences are properly split
        sequences = content.split('>')[1:]  # Skip empty first element
        for seq in sequences:
            header, sequence = seq.strip().split('\n')
            self.assertEqual(len(sequence), 4)  # Each part should be 4 nucleotides
    
    def test_process_fasta_with_mutations(self):
        output_path = os.path.join(self.temp_dir.name, 'output.fasta')
        
        # Test with 50% include and 1.0 mutation rate
        process_fasta(self.test_fasta_path, output_path, 1.0, 50.0)
        
        with open(output_path, 'r') as f:
            content = f.read()
        
        # Check if all sequences are mutated
        sequences = content.split('>')[1:]
        for seq in sequences:
            header, sequence = seq.strip().split('\n')
            original_part = 'ATCG' if 'test_seq1' in header else 'GCTA'
            self.assertNotEqual(sequence, original_part)
            self.assertEqual(len(sequence), 4)

if __name__ == '__main__':
    unittest.main() 