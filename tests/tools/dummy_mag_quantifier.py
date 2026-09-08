#!/usr/bin/env python3
"""Dummy MAG quantifier: one count per MAG FASTA."""

from __future__ import annotations

import argparse
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", required=True)
    parser.add_argument("-r", required=True)
    parser.add_argument("-o", required=True)
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.parent.mkdir(parents=True, exist_ok=True)
    mag_dir = Path(args.r)
    rows = ["mag_id\tN"]
    if mag_dir.is_dir():
        for path in sorted(mag_dir.iterdir()):
            if path.suffix.lower() in {".fa", ".fna", ".fasta"}:
                rows.append(f"{path.stem}\t10")
    elif mag_dir.is_file():
        rows.append(f"{mag_dir.stem}\t10")
    if len(rows) == 1:
        rows.append("mag1\t10")
    dest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
