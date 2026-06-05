#!/usr/bin/env python3

import argparse
import gzip
import os
from collections import Counter
from itertools import product

import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import IncrementalPCA


# ----------------------------
# FASTQ reader
# ----------------------------
def read_fastq(path):
    """Generator to read FASTQ files, handling both raw and gzipped formats."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as f:
        while True:
            header = f.readline().strip()
            if not header:
                break
            seq = f.readline().strip()
            plus = f.readline()
            qual = f.readline().strip()

            # [FIX]: Clean read ID to match the annotation table exactly (remove paired-end suffixes)
            raw_id = header[1:].split()[0]
            clean_id = raw_id.replace("/1", "").replace("/2", "")

            yield clean_id, seq, qual


# ----------------------------
# Entropy metrics
# ----------------------------
def shannon_entropy(seq):
    counts = Counter(seq)
    probs = np.array(list(counts.values())) / len(seq)
    return -np.sum(probs * np.log2(probs + 1e-12))


def renyi_entropy(seq, alpha=2):
    counts = Counter(seq)
    probs = np.array(list(counts.values())) / len(seq)
    return (1 / (1 - alpha)) * np.log2(np.sum(probs**alpha) + 1e-12)


def min_entropy(seq):
    counts = Counter(seq)
    probs = np.array(list(counts.values())) / len(seq)
    return -np.log2(np.max(probs) + 1e-12)


# ----------------------------
# Sequence statistics
# ----------------------------
def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return (g + c) / len(seq)


def gc_skew(seq):
    g = seq.count("G")
    c = seq.count("C")
    return (g - c) / (g + c + 1e-12)


def at_skew(seq):
    a = seq.count("A")
    t = seq.count("T")
    return (a - t) / (a + t + 1e-12)


def mean_quality(qual):
    return np.mean([ord(c) - 33 for c in qual])


# ----------------------------
# k-mer features
# ----------------------------
def kmer_counts(seq, k=4):
    kmers = ["".join(p) for p in product("ACGT", repeat=k)]
    counts = dict.fromkeys(kmers, 0)

    for i in range(len(seq) - k + 1):
        kmer = seq[i : i + k]
        if set(kmer) <= set("ACGT"):
            counts[kmer] += 1

    total = sum(counts.values()) + 1e-12
    return np.array([counts[k] / total for k in kmers])


# ----------------------------
# Complexity metrics
# ----------------------------
def linguistic_complexity(seq, k=4):
    observed = set(seq[i : i + k] for i in range(len(seq) - k + 1))
    possible = min(4**k, len(seq) - k + 1)
    return len(observed) / (possible + 1e-12)


def repeat_fraction(seq, window=10):
    repeats = 0
    total = len(seq) - window + 1
    for i in range(total):
        sub = seq[i : i + window]
        if len(set(sub)) < window / 2:
            repeats += 1
    return repeats / (total + 1e-12)


def dust_score(seq, window=3):
    score = 0
    for i in range(len(seq) - window + 1):
        sub = seq[i : i + window]
        counts = Counter(sub)
        score += sum(v * (v - 1) for v in counts.values())
    return score / (len(seq) + 1e-12)


# ----------------------------
# Main (Two-Pass Batched version)
# ----------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Extract features from FASTQ files safely using Two-Pass chunking."
    )
    parser.add_argument("fastq", help="Input FASTQ file")
    parser.add_argument(
        "-o", "--output", default="features.tsv", help="Output TSV file"
    )
    parser.add_argument("--k", type=int, default=4, help="k-mer size")
    parser.add_argument(
        "--clusters", type=int, default=5, help="Number of KMeans clusters"
    )
    parser.add_argument(
        "--chunk_size",
        type=int,
        default=50000,
        help="Number of reads to process in RAM at once",
    )
    args = parser.parse_args()

    # Initialize online models for dimensionality reduction and clustering
    ipca = IncrementalPCA(n_components=5)
    mbk = MiniBatchKMeans(
        n_clusters=args.clusters,
        random_state=42,
        batch_size=args.chunk_size,
        n_init="auto",
    )

    # Remove existing output file to avoid appending to old data
    if os.path.exists(args.output):
        os.remove(args.output)

    # ---------------------------------------------------------
    # PASS 1: Train Models (IncrementalPCA & MiniBatchKMeans)
    # ---------------------------------------------------------
    print(
        f"[INFO] PASS 1: Training models on {args.fastq} in chunks of {args.chunk_size} reads..."
    )
    kmer_matrix = []
    processed_count = 0

    for _, seq, _ in read_fastq(args.fastq):
        seq = seq.upper()
        kmer_matrix.append(kmer_counts(seq, args.k))
        processed_count += 1

        if len(kmer_matrix) >= args.chunk_size:
            k_mat = np.array(kmer_matrix)
            ipca.partial_fit(k_mat)
            mbk.partial_fit(k_mat)
            kmer_matrix = []
            print(f"  -> Trained on {processed_count} reads...")

    # Tail chunk for Pass 1
    if kmer_matrix:
        k_mat = np.array(kmer_matrix)
        ipca.partial_fit(k_mat)
        mbk.partial_fit(k_mat)
        print(f"  -> Trained on {processed_count} reads (completed Pass 1).")

    # ---------------------------------------------------------
    # PASS 2: Extract Features, Transform, and Write to Disk
    # ---------------------------------------------------------
    print(
        f"\n[INFO] PASS 2: Extracting biological features and writing to {args.output}..."
    )
    records = []
    kmer_matrix = []
    processed_count = 0
    is_first_chunk = True

    for rid, seq, qual in read_fastq(args.fastq):
        seq = seq.upper()

        feats = {
            "read_id": rid,
            "length": len(seq),
            "gc": gc_content(seq),
            "gc_skew": gc_skew(seq),
            "at_skew": at_skew(seq),
            "shannon": shannon_entropy(seq),
            "renyi2": renyi_entropy(seq, 2),
            "min_entropy": min_entropy(seq),
            "ling_complexity": linguistic_complexity(seq, args.k),
            "repeat_frac": repeat_fraction(seq),
            "dust": dust_score(seq),
            "mean_qual": mean_quality(qual),
        }

        kmer_matrix.append(kmer_counts(seq, args.k))
        records.append(feats)
        processed_count += 1

        # Process chunk and write to disk
        if len(records) >= args.chunk_size:
            df = pd.DataFrame(records)
            k_mat = np.array(kmer_matrix)

            # Transform using fully trained models
            kmer_pca = ipca.transform(k_mat)
            for i in range(kmer_pca.shape[1]):
                df[f"kmer_pca_{i+1}"] = kmer_pca[:, i]

            df["tetrabin"] = mbk.predict(k_mat)

            # Write to disk
            df.to_csv(
                args.output, sep="\t", index=False, mode="a", header=is_first_chunk
            )

            # Clear memory
            records = []
            kmer_matrix = []
            is_first_chunk = False
            print(f"  -> Processed and saved {processed_count} reads...")

    # Tail chunk for Pass 2
    if records:
        df = pd.DataFrame(records)
        k_mat = np.array(kmer_matrix)

        kmer_pca = ipca.transform(k_mat)
        for i in range(kmer_pca.shape[1]):
            df[f"kmer_pca_{i+1}"] = kmer_pca[:, i]

        df["tetrabin"] = mbk.predict(k_mat)
        df.to_csv(args.output, sep="\t", index=False, mode="a", header=is_first_chunk)
        print(f"  -> Processed and saved {processed_count} reads (completed Pass 2).")

    print(f"\n[SUCCESS] Features are successfully saved to {args.output}")


if __name__ == "__main__":
    main()
