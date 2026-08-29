#!/usr/bin/env python3

import os
import random
from typing import List, Optional, Tuple
from pathlib import Path
import argparse

from samovar.seqio import (
    fasta_handle,
    is_fasta_name,
    open_text,
    taxid_from_fasta_name,
)


def read_fasta(file_path: str) -> List[Tuple[str, str]]:
    """Read FASTA file (plain or gzipped) and return (header, sequence) tuples."""
    sequences = []
    current_header = None
    current_sequence = []

    with fasta_handle(file_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    sequences.append((current_header, "".join(current_sequence)))
                current_header = line[1:]
                current_sequence = []
            else:
                current_sequence.append(line)

    if current_header is not None:
        sequences.append((current_header, "".join(current_sequence)))

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

def preprocess_fasta(
    input_file: str,
    output_file: str,
    mutation_rate: float,
    include_percent: float,
    inject_taxid: Optional[bool] = None,
):
    """Process FASTA file, apply mutations and split sequences."""
    sequences = read_fasta(input_file)
    input_filename = taxid_from_fasta_name(input_file) or Path(
        str(input_file).removesuffix(".gz")
    ).stem
    if inject_taxid is None:
        flag = os.environ.get("SAMOVAR_INJECT_TAXID", "1").strip().lower()
        inject_taxid = flag not in {"0", "false", "no", "off"}

    out_parent = Path(output_file).parent
    if str(out_parent):
        out_parent.mkdir(parents=True, exist_ok=True)

    with open_text(output_file, "wt") as f:
        for i, (header, sequence) in enumerate(sequences):
            # Calculate length of sequence to include
            include_length = int(len(sequence) * include_percent / 100)
            if include_length <= 0:
                include_length = len(sequence) or 1

            # Split sequence into parts
            for j in range(0, len(sequence), include_length):
                part = sequence[j:j + include_length]
                if not part:  # Skip empty parts
                    continue
                    
                # Apply mutations
                mutated_part = apply_mutations(part, mutation_rate)
                
                # Write to output file
                if inject_taxid:
                    new_header = f">{input_filename}|taxid:{input_filename}|{i+1}|{j//include_length + 1}"
                else:
                    token = header.split()[0] if header else input_filename
                    new_header = f">{token}"
                f.write(f"{new_header}\n")
                f.write(f"{mutated_part}\n")

def process_fasta_directories(directories: List[str]) -> dict:
    """
    Process FASTA files from input directories and extract taxids from filenames.
    Handles various FASTA extensions and gzipped files without decompressing them.
    Also scans ``processed/`` subdirectories (samovar genomes layout).
    """
    result = {}

    def add_file(fasta_file: Path) -> None:
        if not fasta_file.is_file():
            return
        if not is_fasta_name(fasta_file.name, protein=True, nucleotide=True):
            return
        taxid = taxid_from_fasta_name(fasta_file)
        if not taxid:
            return
        from samovar.genome_index import numeric_taxid_for

        resolved = numeric_taxid_for(taxid)
        result[str(fasta_file)] = resolved if resolved else taxid

    for directory in directories:
        dir_path = Path(directory)
        if not dir_path.exists():
            continue
        search = [dir_path]
        nested = dir_path / "processed"
        if nested.is_dir():
            search.append(nested)
        for folder in search:
            try:
                children = list(folder.iterdir())
            except OSError:
                continue
            for fasta_file in children:
                add_file(fasta_file)

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
