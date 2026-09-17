#!/usr/bin/env python3
"""Per-read GC fraction Feature extractor.

Reads FASTQ via the annotator CLI (``-i``/``-I``/``-o``/``-t``) and writes a
headered TSV: ``seq`` plus ``GC`` (ACGT GC fraction; R1 and R2 are concatenated).
This is a Feature sub-contract of Annotation: no taxID column.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Dict, Optional

from samovar.kmer2 import iter_fastq_records
from samovar.qc import gc_fraction


def extract_features(r1: str, r2: Optional[str] = None) -> Dict[str, str]:
    """Map read ID → concatenated sequence (R1 then R2)."""
    by_id: Dict[str, str] = {}
    for path in (r1, r2):
        if not path:
            continue
        for read_id, seq in iter_fastq_records(path):
            if not read_id:
                continue
            if read_id not in by_id:
                by_id[read_id] = seq
            else:
                by_id[read_id] += seq
    return by_id


def write_table(by_id: Dict[str, str], dest: str) -> int:
    path = Path(dest)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as out:
        out.write("seq\tGC\n")
        for read_id, seq in by_id.items():
            out.write(f"{read_id}\t{gc_fraction(seq):.6g}\n")
    return len(by_id)


def parse_output(file_path: str):
    from samovar.parse_annotators import read_custom_raw

    return read_custom_raw(file_path)


def main(argv: Optional[list] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Per-read GC fraction for the Feature contract."
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
