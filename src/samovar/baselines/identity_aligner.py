#!/usr/bin/env python3
"""Built-in aligner baseline: map R1/R2 ids onto MAG FASTA stems under -r."""

from __future__ import annotations

import argparse
import re
import zlib
from pathlib import Path

_TAXID = re.compile(r"taxid[:_](\d+)", re.I)


def _fastq_ids(path: str):
    p = Path(path)
    if not p.is_file():
        return
    with p.open() as handle:
        for i, line in enumerate(handle):
            if i % 4 == 0:
                yield line[1:].split()[0].strip()


def _mag_names(reference: str) -> list[str]:
    src = Path(reference)
    names: list[str] = []
    if src.is_dir():
        for path in sorted(src.iterdir()):
            if path.suffix.lower() in {".fa", ".fna", ".fasta"}:
                names.append(path.stem)
    elif src.is_file() and src.suffix.lower() in {".fa", ".fna", ".fasta"}:
        names.append(src.stem)
    return names or ["mag1"]


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
    mags = _mag_names(args.r)
    lines = ["@HD\tVN:1.6\tSO:unsorted"]
    for mag in mags:
        lines.append(f"@SQ\tSN:{mag}\tLN:4")
    for seq in list(_fastq_ids(args.i)) + list(_fastq_ids(args.r2)):
        tax = _TAXID.search(seq)
        if tax:
            mag = mags[int(tax.group(1)) % len(mags)]
        else:
            mag = mags[zlib.crc32(seq.encode("utf-8")) % len(mags)]
        lines.append(f"{seq}\t0\t{mag}\t1\t60\t4M\t*\t0\t0\tACGT\tIIII")
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
