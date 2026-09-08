#!/usr/bin/env python3
"""Dummy aligner: writes a SAM (named .bam) mapping R1 ids to mag1."""

from __future__ import annotations

import argparse
from pathlib import Path


def _fastq_ids(path: str):
    p = Path(path)
    if not p.is_file():
        return
    with p.open() as handle:
        for i, line in enumerate(handle):
            if i % 4 == 0:
                yield line[1:].split()[0].strip()


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True)
    parser.add_argument("-I", dest="r2", default="")
    parser.add_argument("-r", required=True)
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.parent.mkdir(parents=True, exist_ok=True)
    lines = ["@HD\tVN:1.6\tSO:unsorted", "@SQ\tSN:mag1\tLN:4"]
    for seq in _fastq_ids(args.i):
        lines.append(
            f"{seq}\t0\tmag1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII"
        )
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
