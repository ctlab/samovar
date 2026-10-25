"""Helpers for the C++ annotation sort-merge combiner."""

from __future__ import annotations

import csv
import fcntl
import os
import re
import subprocess
from pathlib import Path
from typing import Optional, Sequence

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


def combine_with_cpp(
    input_dir: str,
    output_dir: str,
    split_n: int,
    chunk_rows: int = 500000,
    truth_table: Optional[str] = None,
) -> None:
    import gzip
    import tempfile

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
    tmp_truth = None
    if truth_table:
        src = Path(truth_table)
        table_path = src
        if is_gzip_file(src):
            tmp_truth = tempfile.NamedTemporaryFile(
                prefix="samovar_truth_", suffix=".tsv", delete=False
            )
            tmp_truth.close()
            with gzip.open(src, "rt", encoding="utf-8", errors="replace") as inn, open(
                tmp_truth.name, "w", encoding="utf-8"
            ) as out:
                out.write(inn.read())
            table_path = Path(tmp_truth.name)
        cmd.extend(["--truth-table", str(table_path)])
    try:
        subprocess.check_call(cmd)
    finally:
        if tmp_truth is not None:
            try:
                Path(tmp_truth.name).unlink()
            except OSError:
                pass


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
    if truth_table:
        cmd.extend(["--truth-table", str(truth_table)])
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
            raise RuntimeError(f"Annotator output has no records: {out_file}")
        if not csv_path.is_file():
            _restore_sample_csv(csv_path, sample)
            if not csv_path.is_file() and _combined_has_sample(annotation_dir, sample):
                raise RuntimeError(
                    f"Cannot add {out_file.name}: {sample} is already in the combined "
                    "table but its columns could not be restored"
                )
        _append_annotator(
            out_file, csv_path, truth_table=truth_table, chunk_rows=chunk_rows
        )
        if not csv_path.is_file() or not _tool_already_ingested(csv_path, out_file.name):
            raise RuntimeError(f"Annotator {out_file.name} was not committed to {csv_path}")
        out_file.unlink()
        sentinel.write_text("ok\n", encoding="utf-8")
    return csv_path


_TOOL_COL = re.compile(r"^(taxID|feat)_(.+)_(\d+)(?:_(.*))?$")
_RESTORE_SKIP = {"seq", "sample", "length", "true", "read_type"}


def _pipe_ids(n: int, text: str) -> list:
    if n <= 0:
        return []
    if text == "" and n == 1:
        return [""]
    parts = text.split("|")
    if len(parts) != n:
        raise RuntimeError("tools sidecar id list length mismatch")
    return parts


def _read_tools(csv_path: Path) -> list:
    side = _sidecar(csv_path)
    if not side.is_file():
        return []
    tools = []
    for line in side.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            raise RuntimeError(f"Bad tools sidecar for {csv_path}")
        n_tax = int(fields[2])
        n_feat = int(fields[3])
        tools.append(
            {
                "filename": fields[0],
                "tool": fields[1],
                "tax_ids": _pipe_ids(n_tax, fields[4] if len(fields) > 4 else ""),
                "feat_ids": _pipe_ids(n_feat, fields[5] if len(fields) > 5 else ""),
            }
        )
    return tools


def _write_tools(csv_path: Path, tools: list) -> None:
    lines = []
    for tool in tools:
        lines.append(
            "\t".join(
                [
                    tool["filename"],
                    tool["tool"],
                    str(len(tool["tax_ids"])),
                    str(len(tool["feat_ids"])),
                    "|".join(tool["tax_ids"]),
                    "|".join(tool["feat_ids"]),
                ]
            )
        )
    text = "\n".join(lines)
    _sidecar(csv_path).write_text(text + ("\n" if text else ""), encoding="utf-8")


def _tools_from_header(header: list) -> list:
    groups = {}
    for name in header:
        if name in _RESTORE_SKIP:
            continue
        match = _TOOL_COL.match(name)
        if not match:
            continue
        prefix, tool, index_s, ident = match.groups()
        key = (int(index_s), tool)
        spec = groups.setdefault(key, {"tax_ids": [], "feat_ids": []})
        bucket = "tax_ids" if prefix == "taxID" else "feat_ids"
        spec[bucket].append("" if ident is None else ident)
    tools = []
    for position, (index_s, tool) in enumerate(sorted(groups)):
        spec = groups[(index_s, tool)]
        tools.append(
            {
                "filename": f"{position:06d}_{tool}.out",
                "tool": tool,
                "tax_ids": spec["tax_ids"],
                "feat_ids": spec["feat_ids"],
            }
        )
    return tools


def _column_name(prefix: str, tool: str, index: int, ident: str) -> str:
    base = f"{prefix}_{tool}_{index}"
    return base if ident == "" else f"{base}_{ident}"


def _combined_has_sample(annotation_dir: Path, sample: str) -> bool:
    combined = annotation_dir / "combined_annotation_table.csv"
    if not combined.is_file():
        return False
    with combined.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            return False
        if "sample" not in reader.fieldnames:
            return True
        return any(row.get("sample") == sample for row in reader)


def _restore_sample_csv(csv_path: Path, sample: str) -> None:
    """Rebuild a deleted per-sample CSV from the combined table before another append.

    Finalize removes the CSV but keeps the tool sidecar. Without this, the next
    annotator starts a new table and the combined file loses the earlier columns.
    """
    if csv_path.is_file():
        return
    combined = csv_path.parent / "combined_annotation_table.csv"
    if not combined.is_file():
        return
    tools = _read_tools(csv_path)
    wrote_sidecar = False
    with combined.open(newline="") as handle:
        reader = csv.DictReader(handle)
        header = list(reader.fieldnames or [])
        if not tools:
            tools = _tools_from_header(header)
            if not tools:
                return
            _write_tools(csv_path, tools)
            wrote_sidecar = True
        columns = ["seq"]
        for index, tool in enumerate(tools):
            for ident in tool["tax_ids"]:
                columns.append(_column_name("taxID", tool["tool"], index, ident))
            for ident in tool["feat_ids"]:
                columns.append(_column_name("feat", tool["tool"], index, ident))
        columns.extend(["length", "true", "read_type"])
        rows = []
        for row in reader:
            if row.get("sample") not in (None, "", sample):
                continue
            rows.append(row)
    if not rows:
        if wrote_sidecar and _sidecar(csv_path).is_file():
            _sidecar(csv_path).unlink()
        return
    with csv_path.open("w", newline="") as out:
        writer = csv.DictWriter(
            out, fieldnames=columns, lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({name: row.get(name, "") for name in columns})


def retrue_annotation(csv_path: Path, truth_table: str) -> None:
    binary = ensure_combine_binary()
    subprocess.check_call(
        [str(binary), "--retrue", str(csv_path), "--truth-table", str(truth_table)]
    )


def finalize_annotations(
    annotation_dir: Path,
    *,
    truth_table: Optional[str] = None,
    delete_parts: bool = True,
) -> Path:
    """Stream per-sample CSVs into ``combined_annotation_table.csv``.

    Per-sample files are removed only after the combined table is in place.
    """
    annotation_dir = Path(annotation_dir)
    parts = sorted(annotation_dir.glob("*.annotation.csv"))
    dest = annotation_dir / "combined_annotation_table.csv"
    if not parts:
        if dest.is_file():
            return dest
        raise FileNotFoundError(f"No per-sample annotation CSVs in {annotation_dir}")
    if truth_table:
        for path in parts:
            retrue_annotation(path, truth_table)
    headers: list[list[str]] = []
    for path in parts:
        with path.open(newline="") as handle:
            header = next(csv.reader(handle), None)
        if not header:
            raise RuntimeError(f"Empty annotation CSV: {path}")
        headers.append(header)
    columns = ["seq"]
    seen = {"seq", "sample"}
    for header in headers:
        for name in header:
            if name in seen or name == "sample":
                continue
            columns.append(name)
            seen.add(name)
    fresh = {path.name[: -len(".annotation.csv")] for path in parts}
    if dest.is_file():
        with dest.open(newline="") as handle:
            previous = next(csv.reader(handle), None) or []
        for name in previous:
            if name in seen or name == "sample":
                continue
            columns.append(name)
            seen.add(name)
    if "sample" not in columns:
        columns.append("sample")
    tmp = dest.with_suffix(".csv.tmp")
    with tmp.open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=columns, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for path in parts:
            sample = path.name[: -len(".annotation.csv")]
            with path.open(newline="") as handle:
                reader = csv.DictReader(handle)
                for row in reader:
                    row["sample"] = sample
                    writer.writerow(row)
        if dest.is_file():
            with dest.open(newline="") as handle:
                reader = csv.DictReader(handle)
                if reader.fieldnames and "sample" in reader.fieldnames:
                    for row in reader:
                        if row.get("sample") in fresh:
                            continue
                        writer.writerow(row)
    tmp.replace(dest)
    if delete_parts:
        for path in parts:
            path.unlink()
            lock = path.with_suffix(".lock")
            if lock.is_file():
                lock.unlink()
    return dest


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
