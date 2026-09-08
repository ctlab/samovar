#!/usr/bin/env python3
"""Dummy MAG taxonomy: MAG dir → mag_id,taxid,lineage."""

from __future__ import annotations

import argparse
from pathlib import Path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", required=True)
    parser.add_argument("-d", dest="db", default="")
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.parent.mkdir(parents=True, exist_ok=True)
    mag_dir = Path(args.c)
    rows = ["mag_id\ttaxid\tlineage"]
    names = []
    if mag_dir.is_dir():
        names = [
            p.stem
            for p in sorted(mag_dir.iterdir())
            if p.suffix.lower() in {".fa", ".fna", ".fasta"}
        ]
    if not names:
        names = ["mag1"]
    for name in names:
        rows.append(f"{name}\t562\td__Bacteria;p__Proteobacteria;g__Escherichia")
    dest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
