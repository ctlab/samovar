#!/usr/bin/env python3
"""Map reads onto an existing MAG directory.

This is the aligner/read-assigner variant that does not invent ids for a new
assembly. A read takes the id of the initial MAG that contains its sequence.
If none does, the built-in aligner assigns a reference name so the table is
still filled.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from samovar.MAG.homology import iter_mag_fastas
from samovar.MAG.tetra import read_fasta_sequence


def _fastq_records(path: str) -> Iterable[Tuple[str, str]]:
    src = Path(path)
    if not path or not src.is_file():
        return
    lines = src.read_text(encoding="utf-8", errors="replace").splitlines()
    for i in range(0, len(lines) - 1, 4):
        if not lines[i].startswith("@"):
            continue
        yield lines[i][1:].split()[0], lines[i + 1].strip().upper()


def _tax_by_mag(path: str) -> Dict[str, str]:
    src = Path(path)
    out: Dict[str, str] = {}
    if not src.is_file():
        return out
    for line in src.read_text(encoding="utf-8").splitlines()[1:]:
        parts = line.split("\t")
        if len(parts) >= 2 and parts[0]:
            out[parts[0]] = parts[1]
    return out


def _mag_sequences(mag_dir: Path) -> List[Tuple[str, str]]:
    loaded = []
    for fasta in iter_mag_fastas(mag_dir):
        loaded.append((fasta.stem, read_fasta_sequence(fasta)))
    return loaded


def contained_mag(read: str, mags: Sequence[Tuple[str, str]]) -> str:
    """MAG id whose sequence contains ``read``, preferring a unique hit."""
    if len(read) < 8:
        return ""
    hits = [name for name, sequence in mags if read in sequence]
    if not hits:
        return ""
    if len(hits) == 1:
        return hits[0]
    lengths = {name: len(sequence) for name, sequence in mags}
    return max(hits, key=lambda name: lengths.get(name, 0))


def _aligner_fallback(r1: str, r2: str, mag_dir: Path, dest: Path) -> List[Tuple[str, str]]:
    from samovar.baselines.identity_aligner import main as align_main

    sam = dest.with_suffix(".sam")
    align_main(["-i", r1, "-I", r2 or "", "-r", str(mag_dir), "-o", str(sam), "-t", "1"])
    rows = []
    for line in sam.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line or line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) >= 3 and parts[2] not in {"*", ""}:
            rows.append((parts[0], parts[2].split("~", 1)[0]))
    return rows


def map_reads(
    r1: str,
    r2: str,
    mag_dir: Path,
    *,
    fallback: bool = True,
    dest_for_fallback: Optional[Path] = None,
) -> List[Tuple[str, str, str]]:
    """Return ``(read_id, sequence, mag_id)`` for reads placed on ``mag_dir``."""
    mags = _mag_sequences(Path(mag_dir))
    placed: List[Tuple[str, str, str]] = []
    missed: List[Tuple[str, str]] = []
    for read_id, sequence in list(_fastq_records(r1)) + list(_fastq_records(r2)):
        mag = contained_mag(sequence, mags)
        if mag:
            placed.append((read_id, sequence, mag))
        else:
            missed.append((read_id, sequence))
    if missed and fallback:
        dest = dest_for_fallback or Path(mag_dir) / "map_fallback.tsv"
        by_id = {read_id: mag for read_id, mag in _aligner_fallback(r1, r2, Path(mag_dir), dest)}
        for read_id, sequence in missed:
            mag = by_id.get(read_id, "")
            if mag:
                placed.append((read_id, sequence, mag))
    return placed


def write_mapping_table(
    r1: str,
    r2: str,
    mag_dir: str,
    taxonomy: str,
    dest: str,
    *,
    feature: str = "MAP_MAG_ID",
) -> Path:
    out = Path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)
    tax = _tax_by_mag(taxonomy)
    rows = [f"seq\ttaxID\t{feature}"]
    for read_id, _sequence, mag in map_reads(r1, r2, Path(mag_dir), dest_for_fallback=out):
        rows.append(f"{read_id}\t{tax.get(mag, '0')}\t{mag}")
    if len(rows) == 1:
        rows.append(f"r0\t0\t")
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return out


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m samovar.MAG.mapping")
    parser.add_argument("-i", required=True)
    parser.add_argument("-I", dest="r2", default="")
    parser.add_argument("-r", required=True, help="Initial MAG directory")
    parser.add_argument("-x", default="", help="MAG taxonomy TSV")
    parser.add_argument("-o", required=True)
    parser.add_argument("-t", dest="threads", default="1")
    args = parser.parse_args(list(argv) if argv is not None else None)
    write_mapping_table(args.i, args.r2, args.r, args.x, args.o)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
