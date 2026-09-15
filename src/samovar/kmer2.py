#!/usr/bin/env python3
"""Dinucleotide (k=2) feature extractor.

Reads FASTQ via the annotator CLI (``-i``/``-I``/``-o``/``-t``) and writes a
headered TSV: ``seq`` plus 16 ACGT 2-mer counts (``AA``, ``AC``, …, ``TT``).
This is a Feature sub-contract of Annotation: no taxID column.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path
from typing import Dict, Iterator, List, Optional, TextIO, Tuple

K = 2
BASES = "ACGT"
KMER_IDS: List[str] = [a + b for a in BASES for b in BASES]
KMER_INDEX = {kmer: i for i, kmer in enumerate(KMER_IDS)}


def _open_text(path: str) -> TextIO:
    gzipped = str(path).endswith(".gz")
    if not gzipped:
        try:
            with open(path, "rb") as raw:
                gzipped = raw.read(2) == b"\x1f\x8b"
        except OSError:
            gzipped = False
    if gzipped:
        return gzip.open(path, "rt")
    return open(path, "r", encoding="utf-8", errors="replace")


def iter_fastq_records(path: str) -> Iterator[Tuple[str, str]]:
    """Yield ``(read_id, sequence)`` from a FASTQ file."""
    p = Path(path)
    if not path or not p.is_file() or p.stat().st_size == 0:
        return
    with _open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline()
            handle.readline()
            handle.readline()
            stripped = header.strip()
            if not stripped:
                continue
            if stripped.startswith("@"):
                stripped = stripped[1:]
            tokens = stripped.split()
            if not tokens:
                continue
            read_id = tokens[0].replace("/1", "").replace("/2", "")
            yield read_id, seq.strip().upper()


def count_kmers(seq: str, k: int = K) -> List[int]:
    """Count ACGT k-mers; skip windows that contain non-ACGT bases."""
    counts = [0] * (4**k)
    if k <= 0 or len(seq) < k:
        return counts
    for i in range(len(seq) - k + 1):
        mer = seq[i : i + k]
        idx = KMER_INDEX.get(mer)
        if idx is not None:
            counts[idx] += 1
    return counts


def extract_features(
    r1: str,
    r2: Optional[str] = None,
) -> Dict[str, List[int]]:
    """Map read ID → 16 dinucleotide counts (R1 and R2 are summed)."""
    by_id: Dict[str, List[int]] = {}
    for path in (r1, r2):
        if not path:
            continue
        for read_id, seq in iter_fastq_records(path):
            if not read_id:
                continue
            add = count_kmers(seq, K)
            if read_id not in by_id:
                by_id[read_id] = add
            else:
                cur = by_id[read_id]
                for i, n in enumerate(add):
                    cur[i] += n
    return by_id


def write_table(by_id: Dict[str, List[int]], dest: str) -> int:
    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        out.write("seq\t" + "\t".join(KMER_IDS) + "\n")
        for read_id in by_id:
            counts = by_id[read_id]
            out.write(read_id + "\t" + "\t".join(str(n) for n in counts) + "\n")
    return len(by_id)


def parse_output(file_path: str):
    """Contract helper: headered seq + k-mer columns."""
    from samovar.parse_annotators import read_custom_raw

    return read_custom_raw(file_path)


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Count ACGT dinucleotides (k=2) per read for the Feature contract."
    )
    parser.add_argument("-i", "-1", dest="r1", required=True, help="R1 FASTQ")
    parser.add_argument("-I", "-2", dest="r2", default="", help="R2 FASTQ (optional)")
    parser.add_argument("-d", dest="db", default="", help="Unused (annotator CLI)")
    parser.add_argument("-o", "--output", dest="o", required=True, help="Output TSV")
    parser.add_argument("-t", dest="threads", default="1", help="Unused (annotator CLI)")
    args = parser.parse_args(argv)
    table = extract_features(args.r1, args.r2 or None)
    write_table(table, args.o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
