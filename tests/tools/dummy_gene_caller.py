#!/usr/bin/env python3
"""Dummy gene caller: contig/MAG FASTA → amino-acid FASTA + GFF."""

from __future__ import annotations

import argparse
from pathlib import Path


def _fastas(src: Path) -> list[Path]:
    if src.is_file():
        return [src]
    if not src.is_dir():
        return []
    out = []
    for path in sorted(src.iterdir()):
        if path.suffix.lower() in {".fa", ".fna", ".fasta"}:
            out.append(path)
    return out


def _translate(seq: str) -> str:
    table = {
        "TTT": "F",
        "TTC": "F",
        "TTA": "L",
        "TTG": "L",
        "TCT": "S",
        "TCC": "S",
        "TCA": "S",
        "TCG": "S",
        "TAT": "Y",
        "TAC": "Y",
        "TAA": "*",
        "TAG": "*",
        "TGT": "C",
        "TGC": "C",
        "TGA": "*",
        "TGG": "W",
        "CTT": "L",
        "CTC": "L",
        "CTA": "L",
        "CTG": "L",
        "CCT": "P",
        "CCC": "P",
        "CCA": "P",
        "CCG": "P",
        "CAT": "H",
        "CAC": "H",
        "CAA": "Q",
        "CAG": "Q",
        "CGT": "R",
        "CGC": "R",
        "CGA": "R",
        "CGG": "R",
        "ATT": "I",
        "ATC": "I",
        "ATA": "I",
        "ATG": "M",
        "ACT": "T",
        "ACC": "T",
        "ACA": "T",
        "ACG": "T",
        "AAT": "N",
        "AAC": "N",
        "AAA": "K",
        "AAG": "K",
        "AGT": "S",
        "AGC": "S",
        "AGA": "R",
        "AGG": "R",
        "GTT": "V",
        "GTC": "V",
        "GTA": "V",
        "GTG": "V",
        "GCT": "A",
        "GCC": "A",
        "GCA": "A",
        "GCG": "A",
        "GAT": "D",
        "GAC": "D",
        "GAA": "E",
        "GAG": "E",
        "GGT": "G",
        "GGC": "G",
        "GGA": "G",
        "GGG": "G",
    }
    aa = []
    dna = "".join(c for c in seq.upper() if c in "ACGT")
    if not dna:
        dna = "ATGAAATAA"
    dna = dna + "A" * ((3 - len(dna) % 3) % 3)
    for i in range(0, len(dna), 3):
        aa.append(table.get(dna[i : i + 3], "X"))
    return "".join(p for p in aa if p != "*") or "M"


def _seq(path: Path) -> str:
    parts = []
    if path.is_file():
        for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
            if not line.startswith(">"):
                parts.append(line.strip())
    return "".join(parts)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("-c", required=True)
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(argv)
    dest = Path(args.o)
    dest.mkdir(parents=True, exist_ok=True)
    srcs = _fastas(Path(args.c)) or [Path(args.c)]
    for fasta in srcs:
        stem = fasta.stem if fasta.name else "genes"
        aa = _translate(_seq(fasta) if fasta.is_file() else "ATGAAATAA")
        (dest / f"{stem}.faa").write_text(f">{stem}_1\n{aa}\n", encoding="utf-8")
        (dest / f"{stem}.gff").write_text(
            f"{stem}\tprodigal\tCDS\t1\t{max(3 * len(aa), 3)}\t.\t+\t0\tID={stem}_1\n",
            encoding="utf-8",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
