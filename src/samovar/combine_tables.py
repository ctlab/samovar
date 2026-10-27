"""Helpers for the C++ annotation sort-merge combiner."""

from __future__ import annotations

import csv
import fcntl
import gzip
import os
import subprocess
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Optional, Sequence

from samovar.parse_annotators import (
    DEFAULT_TAX_RANK,
    ensure_taxid_rank_map,
    canonical_taxid,
    taxid_value_columns,
)
from samovar.seqio import is_gzip_file


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
    from samovar.paths import cpp_binary_loadable

    binary = combine_binary()
    source = repo_root() / "src" / "cpp" / "combine_annotations.cpp"
    makefile = repo_root() / "src" / "cpp" / "Makefile"
    if not source.exists():
        raise FileNotFoundError(f"C++ combiner source not found: {source}")
    if binary.exists() and not cpp_binary_loadable(binary):
        binary.unlink()
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


@contextmanager
def _truth_table_for_cpp(truth_table: Optional[str]) -> Iterator[Optional[str]]:
    """Path the C++ combiner can read. Gzip is unpacked; the combiner does not."""
    if not truth_table:
        yield None
        return
    src = Path(truth_table)
    if not is_gzip_file(src):
        yield str(src)
        return
    tmp = tempfile.NamedTemporaryFile(prefix="samovar_truth_", suffix=".tsv", delete=False)
    tmp.close()
    path = Path(tmp.name)
    try:
        with gzip.open(src, "rt", encoding="utf-8", errors="replace") as inn, path.open(
            "w", encoding="utf-8"
        ) as out:
            while True:
                chunk = inn.read(1 << 20)
                if not chunk:
                    break
                out.write(chunk)
        yield str(path)
    finally:
        path.unlink(missing_ok=True)


def combine_with_cpp(
    input_dir: str,
    output_dir: str,
    split_n: int,
    chunk_rows: int = 500000,
    truth_table: Optional[str] = None,
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
    with _truth_table_for_cpp(truth_table) as table_path:
        if table_path:
            cmd.extend(["--truth-table", table_path])
        subprocess.check_call(cmd)


def sample_name_from(filename: str, split_n: int = 1) -> str:
    """Sample id from ``{sample}_{run}.{tool}.out``, matching the C++ combiner."""
    name = Path(filename).name
    last = name.rfind(".")
    prev = name.rfind(".", 0, last) if last > 0 else -1
    if prev > 0:
        stem = name[:prev]
        us = stem.rfind("_")
        if us > 0:
            return stem[:us]
        if stem:
            return stem
    parts = name.split("_")
    n = max(int(split_n), 1)
    if len(parts) < n:
        return name
    return "_".join(parts[:n])


def annotation_csv_path(annotation_dir: Path, sample: str) -> Path:
    return Path(annotation_dir) / f"{sample}.annotation.csv"


def _sidecar(csv_path: Path) -> Path:
    return Path(str(csv_path) + ".tools")


def _tool_already_ingested(csv_path: Path, filename: str) -> bool:
    side = _sidecar(csv_path)
    if not side.is_file():
        return False
    name = Path(filename).name
    for line in side.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.split("\t", 1)[0] == name:
            return True
    return False


def _has_data_lines(path: Path) -> bool:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if text and not text.startswith("#"):
                return True
    return False


def _append_annotator(
    out_file: Path,
    csv_path: Path,
    *,
    truth_table: Optional[str] = None,
    chunk_rows: int = 500000,
) -> None:
    binary = ensure_combine_binary()
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(binary),
        "--add",
        str(out_file),
        "--name",
        out_file.name,
        "--out",
        str(csv_path),
        "--chunk-rows",
        str(chunk_rows),
    ]
    if csv_path.is_file():
        cmd.extend(["--base", str(csv_path)])
    # Per-sample locks do not protect other samples sharing this directory.
    # The native combiner already accepts a private sort/merge workspace.
    with tempfile.TemporaryDirectory(prefix=".combine_", dir=csv_path.parent) as scratch:
        cmd.extend(["--tmp", scratch])
        with _truth_table_for_cpp(truth_table) as table_path:
            if table_path:
                cmd.extend(["--truth-table", table_path])
            subprocess.check_call(cmd)


def ingest_annotator(
    out_file: Path,
    annotation_dir: Path,
    sentinel: Path,
    *,
    split_n: int = 1,
    truth_table: Optional[str] = None,
    chunk_rows: int = 500000,
) -> Path:
    """Fold one per-read annotator file into the sample CSV, then drop that file.

    The sample CSV is replaced only after the combiner exits 0. A parse failure
    leaves both the CSV and the per-read file in place and does not write the sentinel.
    """
    out_file = Path(out_file)
    annotation_dir = Path(annotation_dir)
    sentinel = Path(sentinel)
    sample = sample_name_from(out_file.name, split_n)
    csv_path = annotation_csv_path(annotation_dir, sample)
    annotation_dir.mkdir(parents=True, exist_ok=True)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    lock_path = annotation_dir / f"{sample}.annotation.lock"
    with lock_path.open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        if not out_file.is_file():
            if _tool_already_ingested(csv_path, out_file.name):
                sentinel.write_text("ok\n", encoding="utf-8")
                return csv_path
            raise FileNotFoundError(f"Annotator output is missing: {out_file}")
        if not _has_data_lines(out_file):
            sentinel.write_text("empty\n", encoding="utf-8")
            return csv_path
        if not csv_path.is_file():
            _restore_sample_csv(csv_path, sample)
        _append_annotator(
            out_file, csv_path, truth_table=truth_table, chunk_rows=chunk_rows
        )
        if not csv_path.is_file() or not _tool_already_ingested(csv_path, out_file.name):
            raise RuntimeError(f"Annotator {out_file.name} was not committed to {csv_path}")
        out_file.unlink()
        sentinel.write_text("ok\n", encoding="utf-8")
    return csv_path


def _restore_sample_csv(csv_path: Path, sample: str) -> None:
    """Ask the combiner to rebuild one sample CSV from the combined table."""
    binary = ensure_combine_binary()
    proc = subprocess.run(
        [str(binary), "--restore", str(csv_path), "--sample", sample],
        check=False,
    )
    if proc.returncode == 0:
        return
    if proc.returncode == 3:
        raise RuntimeError(
            f"Cannot add a tool for {sample}: it is already in the combined "
            "table but its columns could not be restored"
        )
    raise subprocess.CalledProcessError(
        proc.returncode, [str(binary), "--restore", str(csv_path), "--sample", sample]
    )


def finalize_annotations(
    annotation_dir: Path,
    *,
    truth_table: Optional[str] = None,
    delete_parts: bool = True,
) -> Path:
    """Stream per-sample CSVs into ``combined_annotation_table.csv`` via the combiner."""
    annotation_dir = Path(annotation_dir)
    binary = ensure_combine_binary()
    cmd = [str(binary), "--finalize", str(annotation_dir)]
    if not delete_parts:
        cmd.append("--keep-parts")
    with _truth_table_for_cpp(truth_table) as table_path:
        if table_path:
            cmd.extend(["--truth-table", table_path])
        subprocess.check_call(cmd)
    return annotation_dir / "combined_annotation_table.csv"


def main(argv: Optional[Sequence[str]] = None) -> int:
    import argparse

    parser = argparse.ArgumentParser(prog="python -m samovar.combine_tables")
    sub = parser.add_subparsers(dest="cmd", required=True)

    ingest = sub.add_parser("ingest")
    ingest.add_argument("--out", required=True)
    ingest.add_argument("--annotations", required=True)
    ingest.add_argument("--sentinel", required=True)
    ingest.add_argument("--split", type=int, default=1)
    ingest.add_argument("--truth-table", default="")
    ingest.add_argument("--chunk-rows", type=int, default=500000)

    final = sub.add_parser("finalize")
    final.add_argument("--annotations", required=True)
    final.add_argument("--truth-table", default="")
    final.add_argument("--keep-parts", action="store_true")

    args = parser.parse_args(argv)
    if args.cmd == "ingest" and not str(args.annotations).strip():
        parser.error("--annotations must be a directory")
    if args.cmd == "ingest":
        ingest_annotator(
            Path(args.out),
            Path(args.annotations),
            Path(args.sentinel),
            split_n=args.split,
            truth_table=args.truth_table or None,
            chunk_rows=args.chunk_rows,
        )
        return 0
    finalize_annotations(
        Path(args.annotations),
        truth_table=args.truth_table or None,
        delete_parts=not args.keep_parts,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
