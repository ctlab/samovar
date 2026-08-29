"""Ground-truth contracts for combined annotation tables.

Default combiner path still parses ``taxid:`` from read/genome headers.
This module is the CAMI-style alternative: a table of sequence ID → taxid.

Accepted tables (auto-detected):

- CAMI reads mapping: ``#anonymous_read_id  genome_id  tax_id  read_id``
- CAMI GSA mapping: ``@@SEQUENCEID  BINID  TAXID ...``
- Two-column TSV/CSV: ``seq`` / ``SEQUENCEID`` + ``taxid`` / ``tax_id`` / ``true``

``from-genomes`` greps FASTA ``>`` headers and *appends* rows so unprepared
NCBI FASTA can be used without renaming contigs to ``taxid:``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple, Union

from samovar.seqio import fasta_handle, is_fasta_name, taxid_from_fasta_name

PathLike = Union[str, Path]

SEQ_NAMES = {
    "anonymous_read_id",
    "sequenceid",
    "seq",
    "read_id",
    "contig_id",
}
TAX_NAMES = {"tax_id", "taxid", "true"}
PARSE_GENOME = "parse-genome"
GROUND_TRUTH_TABLE = "ground-truth-table"
REGENERATED_MODES = (PARSE_GENOME, GROUND_TRUTH_TABLE)


def normalize_regenerated_mode(value: Optional[str]) -> str:
    raw = str(value or PARSE_GENOME).strip().lower().replace("_", "-")
    aliases = {
        "parse-genome": PARSE_GENOME,
        "parsegenome": PARSE_GENOME,
        "header": PARSE_GENOME,
        "headers": PARSE_GENOME,
        "taxid": PARSE_GENOME,
        GROUND_TRUTH_TABLE: GROUND_TRUTH_TABLE,
        "ground-truth": GROUND_TRUTH_TABLE,
        "cami": GROUND_TRUTH_TABLE,
        "table": GROUND_TRUTH_TABLE,
        "gt": GROUND_TRUTH_TABLE,
    }
    if raw not in aliases:
        raise ValueError(
            f"Unknown regenerated-metagenomes={value!r}. "
            f"Use {PARSE_GENOME} or {GROUND_TRUTH_TABLE}."
        )
    return aliases[raw]


def _open_text(path: PathLike):
    target = Path(path)
    name = target.name.lower()
    if name.endswith(".gz"):
        return gzip.open(target, "rt", encoding="utf-8", errors="replace")
    return target.open("r", encoding="utf-8", errors="replace")


def _norm_col(name: str) -> str:
    text = str(name or "").strip().lstrip("#@").strip()
    return text.lower().replace(" ", "_")


def _split_row(line: str) -> List[str]:
    if "\t" in line:
        return line.rstrip("\n\r").split("\t")
    if "," in line and line.count(",") >= 1:
        return next(csv.reader([line.rstrip("\n\r")]))
    return line.rstrip("\n\r").split()


def iter_truth_pairs(path: PathLike) -> Iterator[Tuple[str, str]]:
    """Yield ``(seq_id, taxid)`` from a CAMI or simple truth table."""
    header: Optional[List[str]] = None
    seq_idx: List[int] = []
    tax_col = -1
    with _open_text(path) as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("@") and not line.startswith("@@"):
                continue
            if header is None:
                looks = (
                    line.startswith("#")
                    or line.startswith("@@")
                    or _norm_col(line.split("\t")[0] if "\t" in line else line.split()[0])
                    in SEQ_NAMES | {"tax_id", "taxid", "true"}
                    or "tax_id" in line.lower()
                    or "anonymous_read" in line.lower()
                    or "sequenceid" in line.lower()
                )
                if looks:
                    fields = _split_row(line.lstrip("#"))
                    if fields and fields[0].startswith("@@"):
                        fields[0] = fields[0][2:]
                    header = [_norm_col(c) for c in fields]
                    for i, col in enumerate(header):
                        if col in SEQ_NAMES and i not in seq_idx:
                            seq_idx.append(i)
                        if col in TAX_NAMES or col.endswith("taxid") or col == "tax_id":
                            tax_col = i
                    if tax_col < 0:
                        for i, col in enumerate(header):
                            if "tax" in col:
                                tax_col = i
                                break
                    if tax_col >= 0 and seq_idx:
                        continue
                header = ["seq", "taxid"]
                seq_idx = [0]
                tax_col = 1
            fields = _split_row(line)
            if tax_col >= len(fields):
                continue
            tax = str(fields[tax_col]).strip()
            if not tax or tax in {"NA", "NaN", "0"}:
                continue
            if not seq_idx:
                yield fields[0].strip(), tax
                continue
            for i in seq_idx:
                if i < len(fields):
                    key = fields[i].strip()
                    if key:
                        yield key, tax


def fasta_header_token(header: str) -> str:
    text = header[1:] if header.startswith(">") else header
    return text.split()[0].strip()


def taxid_from_header(header: str, fallback: str = "") -> str:
    token = fasta_header_token(header)
    lower = token.lower()
    idx = lower.find("taxid:")
    if idx >= 0:
        digits = []
        for ch in token[idx + 6 :]:
            if ch.isdigit():
                digits.append(ch)
            else:
                break
        if digits:
            return "".join(digits)
    return str(fallback or "").strip()


def iter_genome_truth_rows(
    fasta: PathLike, *, taxid: str = ""
) -> Iterator[Tuple[str, str]]:
    path = Path(fasta)
    fallback = taxid or taxid_from_fasta_name(path) or ""
    with fasta_handle(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            token = fasta_header_token(line)
            tid = taxid_from_header(line, fallback)
            if token and tid:
                yield token, tid


def append_truth_rows(dest: PathLike, rows: Iterable[Tuple[str, str]]) -> int:
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    new_file = not dest_path.is_file() or dest_path.stat().st_size == 0
    n = 0
    with dest_path.open("a", encoding="utf-8") as handle:
        if new_file:
            handle.write("seq\ttaxid\n")
        for seq, tax in rows:
            seq = str(seq).strip()
            tax = str(tax).strip()
            if not seq or not tax:
                continue
            handle.write(f"{seq}\t{tax}\n")
            n += 1
    return n


def append_from_fasta(fasta: PathLike, dest: PathLike, *, taxid: str = "") -> int:
    return append_truth_rows(dest, iter_genome_truth_rows(fasta, taxid=taxid))


def append_from_genome_dir(genome_dir: PathLike, dest: PathLike) -> int:
    root = Path(genome_dir)
    if not root.is_dir():
        return 0
    n = 0
    folders = [root]
    for extra in ("processed",):
        sub = root / extra
        if sub.is_dir():
            folders.append(sub)
    seen = set()
    for folder in folders:
        for path in sorted(folder.iterdir()):
            if not path.is_file() or not is_fasta_name(path.name, protein=False):
                continue
            key = str(path.resolve()) if path.exists() else str(path)
            if key in seen:
                continue
            seen.add(key)
            n += append_from_fasta(path, dest)
    return n


def write_normalized_table(src: PathLike, dest: PathLike) -> int:
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with dest_path.open("w", encoding="utf-8") as handle:
        handle.write("seq\ttaxid\n")
        for seq, tax in iter_truth_pairs(src):
            handle.write(f"{seq}\t{tax}\n")
            n += 1
    return n


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m samovar.ground_truth")
    sub = parser.add_subparsers(dest="command", required=True)
    genomes = sub.add_parser(
        "from-genomes",
        help="Grep FASTA '>' headers and append seq\\ttaxid rows",
    )
    genomes.add_argument("--genome-dir", "--genome_dir", dest="genome_dir", required=True)
    genomes.add_argument("--output", "-o", required=True)
    genomes.add_argument(
        "--append",
        action="store_true",
        default=True,
        help="Append to the output table (default)",
    )
    genomes.add_argument(
        "--overwrite",
        action="store_true",
        help="Truncate the output table before appending genome headers",
    )
    fasta = sub.add_parser("from-fasta", help="Append one FASTA's headers to a truth table")
    fasta.add_argument("fasta")
    fasta.add_argument("--output", "-o", required=True)
    fasta.add_argument("--taxid", default="")
    norm = sub.add_parser("normalize", help="Rewrite a CAMI/custom table to seq\\ttaxid")
    norm.add_argument("src")
    norm.add_argument("--output", "-o", required=True)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parser().parse_args(list(argv) if argv is not None else None)
    if args.command == "from-genomes":
        dest = Path(args.output)
        if getattr(args, "overwrite", False) and dest.exists():
            dest.unlink()
        n = append_from_genome_dir(args.genome_dir, args.output)
        print(f"Appended {n} header(s) from {args.genome_dir} -> {args.output}")
        return 0
    if args.command == "from-fasta":
        n = append_from_fasta(args.fasta, args.output, taxid=args.taxid)
        print(f"Appended {n} header(s) from {args.fasta} -> {args.output}")
        return 0
    if args.command == "normalize":
        n = write_normalized_table(args.src, args.output)
        print(f"Wrote {n} row(s) to {args.output}")
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
