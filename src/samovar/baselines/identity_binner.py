#!/usr/bin/env python3
"""Built-in binner baseline: split the assembly into two MAG FASTA files."""

from __future__ import annotations

import argparse
from pathlib import Path


def _assembly_seq(src: Path) -> str:
    if not src.is_file():
        return "ACGT"
    seqs = []
    for line in src.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            continue
        seqs.append(line.strip())
    return "".join(seqs) or "ACGT"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", required=True)
    parser.add_argument("-i", dest="r1", default="")
    parser.add_argument("-I", dest="r2", default="")
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.mkdir(parents=True, exist_ok=True)
    seq = _assembly_seq(Path(args.c))
    if len(seq) >= 8:
        mid = len(seq) // 2
        parts = [("mag1", seq[:mid]), ("mag2", seq[mid:])]
    else:
        parts = [("mag1", seq), ("mag2", seq)]
    for name, chunk in parts:
        (dest / f"{name}.fa").write_text(f">{name}\n{chunk}\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
