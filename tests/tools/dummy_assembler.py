#!/usr/bin/env python3
"""Dummy assembler: FASTQ → one contig FASTA."""

from __future__ import annotations

import argparse
from pathlib import Path


def _seq_from_fastq(path: str) -> str:
    p = Path(path)
    if not p.is_file() or p.stat().st_size == 0:
        return "ACGT"
    lines = p.read_text().splitlines()
    seqs = [lines[i] for i in range(1, len(lines), 4) if i < len(lines)]
    return "".join(seqs) or "ACGT"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-I", dest="r2", default="")
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.parent.mkdir(parents=True, exist_ok=True)
    seq = _seq_from_fastq(args.i)
    dest.write_text(f">contig1\n{seq}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
