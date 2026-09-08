#!/usr/bin/env python3
"""Dummy binner: contigs FASTA → MAG directory with one MAG."""

from __future__ import annotations

import argparse
from pathlib import Path


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
    src = Path(args.c)
    text = src.read_text(encoding="utf-8") if src.is_file() else ">contig1\nACGT\n"
    (dest / "mag1.fa").write_text(text, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
