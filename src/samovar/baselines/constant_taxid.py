#!/usr/bin/env python3
"""Built-in annotator: assign one NCBI taxID to every FASTQ read.

Default taxID is 9606 (Homo sapiens). Output is two tab-separated columns
(seq, taxID) in the custom-annotator layout. ``parse_output`` reads the same
shape for ``samovar tools import --type annotator``.
"""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path
from typing import Iterator, Optional, TextIO

import pandas as pd


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
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return
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


def parse_output(file_path: str) -> pd.DataFrame:
    path = Path(file_path)
    if not path.is_file() or path.stat().st_size == 0:
        return pd.DataFrame(columns=["seq", "taxID"])
    df = pd.read_table(path, header=None)
    df = df.iloc[:, :2]
    df.columns = ["seq", "taxID"]
    return df


def classify_fastq(r1: str, output: str, taxid: str = "9606", r2: Optional[str] = None) -> int:
    """Write seq\\ttaxID for every unique read ID in r1 (and r2 if given)."""
    seen = set()
    n = 0
    dest = Path(output)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as out:
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


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Assign a constant NCBI taxID to every FASTQ read."
    )
    parser.add_argument("-i", "-1", dest="r1", required=True, help="R1 FASTQ")
    parser.add_argument("-I", "-2", dest="r2", default=None, help="R2 FASTQ (optional)")
    parser.add_argument("-d", dest="db", default="", help="Unused; accepted for annotator CLI")
    parser.add_argument("-o", "--output", required=True, help="Output TSV (seq, taxID)")
    parser.add_argument("-t", dest="threads", default="1")
    parser.add_argument("--taxid", default="9606", help="NCBI taxID to assign (default: 9606)")
    args = parser.parse_args(argv)
    classify_fastq(args.r1, args.output, taxid=str(args.taxid), r2=args.r2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
