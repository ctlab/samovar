#!/usr/bin/env python3

import random
from typing import List, Tuple
from pathlib import Path
import argparse

def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    """Read FASTA file and return list of (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_sequence = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            if line.startswith('>'):
                if current_header is not None:
                    sequences.append((current_header, ''.join(current_sequence)))
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)
    
    if current_header is not None:
        sequences.append((current_header, ''.join(current_sequence)))
    
    return sequences

def apply_mutations(sequence: str, mutation_rate: float) -> str:
    """Apply random mutations to sequence based on mutation rate."""
    nucleotides = ['A', 'T', 'G', 'C']
    mutated_sequence = list(sequence)
    
    # If mutation rate is 1.0, ensure at least one mutation occurs
    if mutation_rate == 1.0:
        # Choose a random position to mutate
        pos = random.randrange(len(mutated_sequence))
        other_nucleotides = [n for n in nucleotides if n != mutated_sequence[pos]]
        mutated_sequence[pos] = random.choice(other_nucleotides)
    else:
        # Normal mutation process
        for i in range(len(mutated_sequence)):
            if random.random() < mutation_rate:
                # Get other nucleotides excluding current one
                other_nucleotides = [n for n in nucleotides if n != mutated_sequence[i]]
                mutated_sequence[i] = random.choice(other_nucleotides)
    
    return ''.join(mutated_sequence)

def preprocess_fasta(input_file: str, output_file: str, mutation_rate: float, include_percent: float):
    """Process FASTA file, apply mutations and split sequences."""
    sequences = read_fasta(input_file)
    input_filename = Path(input_file).stem
    
    with open(output_file, 'w') as f:
        for i, (header, sequence) in enumerate(sequences):
            # Calculate length of sequence to include
            include_length = int(len(sequence) * include_percent / 100)
            
            # Split sequence into parts
            for j in range(0, len(sequence), include_length):
                part = sequence[j:j + include_length]
                if not part:  # Skip empty parts
                    continue
                    
                # Apply mutations
                mutated_part = apply_mutations(part, mutation_rate)
                
                # Write to output file
                new_header = f">{input_filename}|taxid:{input_filename}|{i+1}|{j//include_length + 1}"
                f.write(f"{new_header}\n")
                f.write(f"{mutated_part}\n")

def process_fasta_directories(directories: List[str]) -> dict:
    """
    Process FASTA files from input directories and extract taxids from filenames.
    Handles various FASTA extensions and gzipped files.
    
    Args:
        directories: List of directory paths containing FASTA files
        
    Returns:
        Dictionary mapping input files to their taxids
    """
    result = {}
    
    for directory in directories:
        dir_path = Path(directory)
        if not dir_path.exists():
            continue
            
        # Look for all possible FASTA files
        for fasta_file in dir_path.glob("*.*"):
            # Handle gzipped files first
            if fasta_file.suffix.lower() == '.gz':
                # Decompress the file
                import subprocess
                try:
                    subprocess.run(['gzip', '-d', str(fasta_file)], check=True)
                    # Update the file path to point to the decompressed file
                    fasta_file = fasta_file.with_suffix('')
                except subprocess.CalledProcessError:
                    continue
            
            # Skip if not a FASTA file
            if not fasta_file.suffix.lower() in ['.fa', '.fna', '.fasta', '.faa', '.frn']:
                continue
            
            # Extract taxid from filename (assuming format like '12345.fa')
            taxid = fasta_file.stem
            if taxid.isdigit():
                result[str(fasta_file)] = taxid
    
    return result

def main():
    parser = argparse.ArgumentParser(description='Process FASTA files with mutations and splitting')
    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('output_file', help='Output FASTA file')
    parser.add_argument('--mutation-rate', type=float, default=0.01,
                      help='Mutation rate (default: 0.01)')
    parser.add_argument('--include-percent', type=float, default=50.0,
                      help='Percentage of sequence to include in each split (default: 50.0)')
    
    args = parser.parse_args()
    preprocess_fasta(args.input_file, args.output_file, args.mutation_rate, args.include_percent)

if __name__ == '__main__':
    main() 