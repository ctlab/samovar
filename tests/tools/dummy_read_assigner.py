#!/usr/bin/env python3
"""Dummy read assigner: SAM/BAM-like text + MAG taxonomy → seq/taxID."""

from __future__ import annotations

import argparse
from pathlib import Path


def _tax_by_mag(path: Path) -> dict:
    out = {}
    if not path.is_file():
        return {"mag1": "562"}
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 2:
            out[parts[0]] = parts[1]
    return out or {"mag1": "562"}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", required=True)
    parser.add_argument("-x", required=True)
    parser.add_argument("-o", required=True)
    args = parser.parse_args(argv)
    tax = _tax_by_mag(Path(args.x))
    dest = Path(args.o)
    dest.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    bam = Path(args.b)
    if bam.is_file():
        for line in bam.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line or line.startswith("@"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            seq, rname = parts[0], parts[2]
            if rname in { "*", "" }:
                continue
            taxid = tax.get(rname, tax.get("mag1", "0"))
            rows.append(f"{seq}\t{taxid}")
    if not rows:
        rows.append("r0\t562")
    dest.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
