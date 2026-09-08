#!/usr/bin/env python3
"""Dummy combine: copy MAG FASTA from the first input directory."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", required=True)
    parser.add_argument("-q", dest="qc", default="")
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.mkdir(parents=True, exist_ok=True)
    first = Path(str(args.c).split(",")[0].strip())
    copied = False
    if first.is_dir():
        for path in sorted(first.iterdir()):
            if path.suffix.lower() in {".fa", ".fna", ".fasta"}:
                shutil.copy2(path, dest / path.name)
                copied = True
    if not copied:
        (dest / "mag1.fa").write_text(">mag1\nACGT\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
