"""Helpers for the C++ annotation sort-merge combiner."""

from __future__ import annotations

import csv
import os
import subprocess
from pathlib import Path

from samovar.parse_annotators import (
    DEFAULT_TAX_RANK,
    ensure_taxid_rank_map,
    canonical_taxid,
    taxid_value_columns,
)


def repo_root() -> Path:
    from samovar.paths import repo_root as _root

    return _root()


def combine_binary() -> Path:
    return repo_root() / "bin" / "samovar_combine_annotations"


def _cxx() -> str:
    from samovar.paths import cxx_compiler

    found = cxx_compiler()
    if not found:
        raise FileNotFoundError(
            "No C++ compiler found (g++, c++, or clang++). "
            "Install a compiler or set CXX."
        )
    return found


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
            env = os.environ.copy()
            env["CXX"] = _cxx()
            subprocess.check_call(["make", "-C", str(makefile.parent)], env=env)
        else:
            binary.parent.mkdir(parents=True, exist_ok=True)
            subprocess.check_call(
                [_cxx(), "-O3", "-std=c++17", "-o", str(binary), str(source)]
            )
    if not os.access(binary, os.X_OK):
        raise RuntimeError(f"C++ combiner is not executable: {binary}")
    return binary


def apply_rank_level_csv(path: str, level: str = DEFAULT_TAX_RANK, cache_path: str = None) -> None:
    """Rewrite a CSV copy's taxID columns using the cached rank map.

    Do not use this on combined annotation tables used for genome fetch.
    """
    src = Path(path)
    with src.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return
        fieldnames = list(reader.fieldnames)
        tax_cols = taxid_value_columns(fieldnames)
        unique = set()
        for row in reader:
            for col in tax_cols:
                value = (row.get(col) or "").strip()
                if value:
                    unique.add(value)

    taxid_map = ensure_taxid_rank_map(unique, level, cache_path)

    tmp = src.with_suffix(src.suffix + ".tmp")
    with src.open(newline="") as fin, tmp.open("w", newline="") as fout:
        reader = csv.DictReader(fin)
        writer = csv.DictWriter(fout, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in reader:
            for col in tax_cols:
                value = row.get(col) or "0"
                key = canonical_taxid(value)
                row[col] = taxid_map.get(value, taxid_map.get(key, key))
            writer.writerow(row)
    tmp.replace(src)


def apply_species_level_csv(path: str, level: str = DEFAULT_TAX_RANK) -> None:
    """Backward-compatible alias for :func:`apply_rank_level_csv`."""
    apply_rank_level_csv(path, level=level)


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
