"""Helpers for the C++ annotation sort-merge combiner."""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

from samovar.parse_annotators import _resolve_taxid_by_rank


def repo_root() -> Path:
    return Path(__file__).resolve().parent.parent.parent


def combine_binary() -> Path:
    return repo_root() / "bin" / "samovar_combine_annotations"


def ensure_combine_binary() -> Path:
    """Build the C++ combiner if the binary is missing or older than the source."""
    binary = combine_binary()
    source = repo_root() / "src" / "cpp" / "combine_annotations.cpp"
    makefile = repo_root() / "src" / "cpp" / "Makefile"
    if not source.exists():
        raise FileNotFoundError(f"C++ combiner source not found: {source}")
    needs_build = (not binary.exists()) or (
        source.stat().st_mtime > binary.stat().st_mtime
    )
    if needs_build:
        if makefile.exists():
            subprocess.check_call(["make", "-C", str(makefile.parent)])
        else:
            binary.parent.mkdir(parents=True, exist_ok=True)
            subprocess.check_call(
                ["g++", "-O3", "-std=c++17", "-o", str(binary), str(source)]
            )
    if not os.access(binary, os.X_OK):
        raise RuntimeError(f"C++ combiner is not executable: {binary}")
    return binary


def apply_species_level_csv(path: str, level: str = "species") -> None:
    """Rewrite taxID_* columns to the requested rank without loading the table."""
    src = Path(path)
    with src.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return
        fieldnames = list(reader.fieldnames)
        tax_cols = [c for c in fieldnames if c.startswith("taxID")]
        unique = set()
        for row in reader:
            for col in tax_cols:
                value = (row.get(col) or "").strip()
                if value:
                    unique.add(value)

    taxid_map = {}
    for taxid in unique:
        if taxid in {"0", "", "nan", "None"}:
            taxid_map[taxid] = "0" if taxid != "" else ""
            continue
        ranked = _resolve_taxid_by_rank(taxid, level)
        taxid_map[taxid] = ranked if ranked is not None else taxid

    tmp = src.with_suffix(src.suffix + ".tmp")
    with src.open(newline="") as fin, tmp.open("w", newline="") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            for col in tax_cols:
                value = row.get(col) or "0"
                row[col] = taxid_map.get(value, value)
            writer.writerow(row)
    tmp.replace(src)


def combine_with_cpp(
    input_dir: str, output_dir: str, split_n: int, chunk_rows: int = 500000
) -> None:
    binary = ensure_combine_binary()
    os.makedirs(output_dir, exist_ok=True)
    cmd = [
        str(binary),
        "-i",
        input_dir,
        "-o",
        output_dir,
        "-s",
        str(split_n),
        "--chunk-rows",
        str(chunk_rows),
    ]
    subprocess.check_call(cmd)
