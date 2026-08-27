#!/usr/bin/env python3
"""Assign a constant NCBI taxID to every FASTQ read (default: Homo sapiens, 9606).

This is a dummy classifier used to test custom annotator integration. Output is
two tab-separated columns: sequence ID and taxID (CustomAnnotator format).
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Iterator, Optional, TextIO


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


def iter_fastq_ids(path: str) -> Iterator[str]:
    """Yield cleaned FASTQ read IDs (paired-end /1 /2 suffixes stripped)."""
    with _open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                break
            handle.readline()
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
            yield tokens[0].replace("/1", "").replace("/2", "")


def classify_fastq(r1: str, output: str, taxid: str = "9606", r2: Optional[str] = None) -> int:
    """Write seq\\ttaxID for every unique read ID in r1 (and r2 if given)."""
    seen = set()
    n = 0
    with open(output, "w", encoding="utf-8") as out:
        for path in (r1, r2):
            if not path or not Path(path).is_file():
                continue
            for read_id in iter_fastq_ids(path):
                if not read_id or read_id in seen:
                    continue
                seen.add(read_id)
                out.write(f"{read_id}\t{taxid}\n")
                n += 1
    return n


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Dummy classifier: assign a constant taxID to all reads."
    )
    parser.add_argument("-i", "-1", dest="r1", required=True, help="R1 FASTQ")
    parser.add_argument("-I", "-2", dest="r2", default=None, help="R2 FASTQ (optional)")
    parser.add_argument("-o", "--output", required=True, help="Output TSV (seq, taxID)")
    parser.add_argument("--taxid", default="9606", help="NCBI taxID to assign (default: 9606)")
    args = parser.parse_args()
    classify_fastq(args.r1, args.output, taxid=str(args.taxid), r2=args.r2)


if __name__ == "__main__":
    main()
