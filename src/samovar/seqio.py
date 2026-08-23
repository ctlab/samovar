"""Flexible FASTA/FASTQ I/O: uncompressed or gzip, several extensions.

Hard-coded ``.fa`` / ``.fq`` suffixes live here so genome prep, database
builders, ISS, and annotator workflows can share the same lookup rules.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import tempfile
from contextlib import contextmanager
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, TextIO, Tuple, Union

PathLike = Union[str, os.PathLike]

NUCLEOTIDE_FASTA_EXTS = (".fasta", ".fna", ".fa", ".fsa", ".ffn")
PROTEIN_FASTA_EXTS = (".faa", ".faa.gz")
GFF_EXTS = (".gff", ".gff3", ".gtf")
FASTQ_EXTS = (".fastq", ".fq")
COMPRESSION_EXTS = (".gz",)

# ISS writes ``prefix_R1.fastq``; we accept those plus .fq and gzip variants.
FASTQ_R1_SUFFIXES = (
    "_R1.fastq.gz",
    "_R1.fq.gz",
    "_R1.fastq",
    "_R1.fq",
)
FASTQ_R2_SUFFIXES = (
    "_R2.fastq.gz",
    "_R2.fq.gz",
    "_R2.fastq",
    "_R2.fq",
)

PROCESSED_FASTA_EXTS = (
    "-processed.fasta.gz",
    "-processed.fna.gz",
    "-processed.fa.gz",
    "-processed.fasta",
    "-processed.fna",
    "-processed.fa",
)


def as_path(path: PathLike) -> Path:
    return Path(path)


def is_gzip_path(path: PathLike) -> bool:
    name = as_path(path).name.lower()
    return name.endswith(".gz") or name.endswith(".bgz") or name.endswith(".gzip")


def strip_compression_suffix(name: str) -> str:
    lower = name.lower()
    for ext in (".gz", ".bgz", ".gzip", ".bz2", ".xz"):
        if lower.endswith(ext):
            return name[: -len(ext)]
    return name


def sequence_stem(path: PathLike) -> str:
    """Filename without compression or sequence-format suffix.

    ``123.fna.gz`` → ``123``; ``562-processed.fasta.gz`` → ``562-processed``.
    """
    name = strip_compression_suffix(as_path(path).name)
    lower = name.lower()
    for ext in NUCLEOTIDE_FASTA_EXTS + (".faa", ".frn") + FASTQ_EXTS + GFF_EXTS:
        if lower.endswith(ext):
            return name[: -len(ext)]
    return Path(name).stem


def taxid_from_fasta_name(path: PathLike) -> Optional[str]:
    """Numeric taxid from ``{taxid}.fa``, ``{taxid}-processed.fasta.gz``, etc."""
    name = strip_compression_suffix(as_path(path).name)
    lower = name.lower()
    for marker in ("-processed.fasta", "-processed.fna", "-processed.fa"):
        if lower.endswith(marker):
            stem = name[: -len(marker)]
            return stem if stem.isdigit() else None
    stem = sequence_stem(path)
    return stem if stem.isdigit() else None


def is_fasta_name(name: str, *, protein: bool = True, nucleotide: bool = True) -> bool:
    base = strip_compression_suffix(name).lower()
    if nucleotide and any(base.endswith(ext) for ext in NUCLEOTIDE_FASTA_EXTS + (".frn",)):
        return True
    if protein and base.endswith(".faa"):
        return True
    return False


def is_fastq_name(name: str) -> bool:
    base = strip_compression_suffix(name).lower()
    return any(base.endswith(ext) for ext in FASTQ_EXTS)


def is_gff_name(name: str) -> bool:
    base = strip_compression_suffix(name).lower()
    return any(base.endswith(ext) for ext in GFF_EXTS)


def open_text(path: PathLike, mode: str = "rt") -> TextIO:
    """Open a text file, using gzip when the path (or content) is compressed."""
    path = as_path(path)
    if "b" in mode:
        raise ValueError("open_text is for text modes only")
    if is_gzip_path(path):
        return gzip.open(path, mode, encoding="utf-8", errors="replace")
    return open(path, mode, encoding="utf-8", errors="replace")


@contextmanager
def fasta_handle(path: PathLike, mode: str = "rt") -> Iterator[TextIO]:
    handle = open_text(path, mode)
    try:
        yield handle
    finally:
        handle.close()


def iter_seqio_fasta(path: PathLike):
    """Yield Biopython SeqRecords from an uncompressed or gzipped FASTA."""
    from Bio import SeqIO

    with fasta_handle(path) as handle:
        yield from SeqIO.parse(handle, "fasta")


def gzip_file(path: PathLike, *, remove_source: bool = True) -> Path:
    """Compress ``path`` to ``path.gz`` (no-op if already gzipped)."""
    src = as_path(path)
    if not src.exists():
        raise FileNotFoundError(src)
    if is_gzip_path(src):
        return src
    dest = Path(str(src) + ".gz")
    with open(src, "rb") as fin, gzip.open(dest, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    if remove_source:
        try:
            src.unlink()
        except OSError:
            pass
    return dest


def gunzip_file(path: PathLike, dest: Optional[PathLike] = None, *, remove_source: bool = False) -> Path:
    src = as_path(path)
    if not is_gzip_path(src):
        return src
    out = as_path(dest) if dest is not None else Path(strip_compression_suffix(str(src)))
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(src, "rb") as fin, open(out, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    if remove_source:
        try:
            src.unlink()
        except OSError:
            pass
    return out


def maybe_gzip(path: PathLike, enabled: bool = True) -> Path:
    src = as_path(path)
    if not enabled or not src.exists():
        return src
    return gzip_file(src)


def processed_genome_path(genome_dir: PathLike, taxid: str, gzip_genomes: bool = True) -> Path:
    taxid = str(taxid).split(".")[0]
    suffix = "-processed.fasta.gz" if gzip_genomes else "-processed.fasta"
    return as_path(genome_dir) / f"{taxid}{suffix}"


def find_existing_processed_genome(genome_dir: PathLike, taxid: str) -> Optional[Path]:
    taxid = str(taxid).split(".")[0]
    directory = as_path(genome_dir)
    for ext in PROCESSED_FASTA_EXTS:
        candidate = directory / f"{taxid}{ext}"
        if candidate.exists():
            return candidate
    return None


def genome_lookup_extensions() -> Tuple[str, ...]:
    """Prefer processed genomes, then gzipped raw, then uncompressed raw."""
    return PROCESSED_FASTA_EXTS + (
        ".fa.gz",
        ".fna.gz",
        ".fasta.gz",
        ".fsa.gz",
        ".fa",
        ".fna",
        ".fasta",
        ".fsa",
    )


def find_genome_file(genome_dir: PathLike, taxid: str) -> Optional[str]:
    taxid = str(taxid).split(".")[0]
    directory = as_path(genome_dir)
    for ext in genome_lookup_extensions():
        candidate = directory / f"{taxid}{ext}"
        if candidate.exists():
            return str(candidate)
    return None


def list_fasta_files(
    directory: PathLike,
    *,
    nucleotide: bool = True,
    protein: bool = False,
) -> List[Path]:
    directory = as_path(directory)
    if not directory.is_dir():
        return []
    found: List[Path] = []
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        if is_fasta_name(path.name, protein=protein, nucleotide=nucleotide):
            found.append(path)
    return found


def mate_suffixes(mate: str) -> Tuple[str, ...]:
    mate = mate.upper()
    if mate == "R1":
        return FASTQ_R1_SUFFIXES
    if mate == "R2":
        return FASTQ_R2_SUFFIXES
    raise ValueError(f"mate must be R1 or R2, not {mate!r}")


def sample_name_from_r1(path: PathLike) -> str:
    name = as_path(path).name
    for suffix in FASTQ_R1_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"Not a recognised R1 FASTQ name: {path}")


def r2_for_r1(r1: PathLike) -> Path:
    path = as_path(r1)
    name = path.name
    for s1, s2 in zip(FASTQ_R1_SUFFIXES, FASTQ_R2_SUFFIXES):
        if name.endswith(s1):
            return path.with_name(name[: -len(s1)] + s2)
    raise ValueError(f"Not a recognised R1 FASTQ name: {r1}")


def find_fastq_mate(directory: PathLike, sample: str, mate: str = "R1") -> Optional[Path]:
    directory = as_path(directory)
    for suffix in mate_suffixes(mate):
        candidate = directory / f"{sample}{suffix}"
        if candidate.exists():
            return candidate
    return None


# CAMISIM harvest writes ``{n}_{tech}_R1.fastq`` then copies to ``{n}_full_``.
# Annotators must only see the ISS-style ``*_full`` names.
_CAMISIM_TECH_SAMPLE = re.compile(
    r"^(?P<n>\d+)_(illumina|ont|wgsim|nanosim3|art)$",
    re.IGNORECASE,
)
_FULL_SAMPLE = re.compile(r"^(?P<n>\d+)_full$", re.IGNORECASE)


def _drop_duplicate_camisim_tech_samples(samples: Sequence[str]) -> List[str]:
    full_ns = {m.group("n") for s in samples if (m := _FULL_SAMPLE.match(s))}
    if not full_ns:
        return list(samples)
    return [
        s
        for s in samples
        if not ((m := _CAMISIM_TECH_SAMPLE.match(s)) and m.group("n") in full_ns)
    ]


def list_fastq_samples(directory: PathLike) -> List[str]:
    directory = as_path(directory)
    if not directory.is_dir():
        return []
    samples: List[str] = []
    seen = set()
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        name = path.name
        for suffix in FASTQ_R1_SUFFIXES:
            if name.endswith(suffix):
                sample = name[: -len(suffix)]
                if sample not in seen:
                    seen.add(sample)
                    samples.append(sample)
                break
    return _drop_duplicate_camisim_tech_samples(samples)


def list_r1_files(directory: PathLike) -> List[Path]:
    directory = as_path(directory)
    files: List[Path] = []
    for sample in list_fastq_samples(directory):
        found = find_fastq_mate(directory, sample, "R1")
        if found is not None:
            files.append(found)
    return files


def fastq_pair_paths(prefix: PathLike, gzip_reads: bool = False) -> Tuple[str, str]:
    """Return ``{prefix}_R1.fastq[.gz]`` and R2 paths (not necessarily existing)."""
    prefix_s = str(prefix)
    ext = ".fastq.gz" if gzip_reads else ".fastq"
    return f"{prefix_s}_R1{ext}", f"{prefix_s}_R2{ext}"


def iter_fastq_records(path: PathLike) -> Iterator[Tuple[str, str, str, str]]:
    with open_text(path) as handle:
        while True:
            header = handle.readline()
            if not header:
                return
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not (seq and plus and qual):
                return
            yield header, seq, plus, qual


def write_text_lines(path: PathLike, chunks: Iterable[str]) -> None:
    dest = as_path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open_text(dest, "wt") as handle:
        wrote = False
        for chunk in chunks:
            handle.write(chunk)
            wrote = True
        if not wrote:
            handle.write("\n")


def concat_fastq_files(sources: Sequence[PathLike], dest: PathLike) -> Path:
    dest_p = as_path(dest)
    dest_p.parent.mkdir(parents=True, exist_ok=True)
    with open_text(dest_p, "wt") as out:
        wrote = False
        for src in sources:
            src_p = as_path(src)
            if not src_p.exists():
                continue
            with open_text(src_p) as handle:
                text = handle.read()
            if text.strip():
                out.write(text)
                if not text.endswith("\n"):
                    out.write("\n")
                wrote = True
        if not wrote:
            out.write("\n")
    return dest_p


def concat_r1_fastqs(directory: PathLike, dest: PathLike) -> Path:
    return concat_fastq_files(list_r1_files(directory), dest)


def link_or_copy_reads(src_dir: PathLike, dst_dir: PathLike) -> int:
    """Symlink (or copy) FASTQ pairs from ``src_dir`` into ``dst_dir``."""
    src = as_path(src_dir).resolve()
    dst = as_path(dst_dir).resolve()
    dst.mkdir(parents=True, exist_ok=True)
    if src == dst:
        return 0
    n = 0
    for r1 in list_r1_files(src):
        r2 = r2_for_r1(r1)
        for src_file in (r1, r2):
            if not src_file.exists():
                continue
            target = dst / src_file.name
            if target.exists() or target.is_symlink():
                continue
            try:
                target.symlink_to(src_file.resolve())
            except OSError:
                shutil.copy2(src_file, target)
            n += 1
    return n


def has_r1_reads(directory: PathLike) -> bool:
    return bool(list_fastq_samples(directory))


@contextmanager
def uncompressed_fasta_for_tool(path: PathLike) -> Iterator[str]:
    """Yield a plain FASTA path; decompress to a temp file only when needed.

    ISS (and a few other CLIs) do not read ``.gz`` genomes. The original
    gzipped file is left in place.
    """
    src = as_path(path)
    if not is_gzip_path(src):
        yield str(src)
        return
    fd, tmp = tempfile.mkstemp(suffix=".fasta", prefix="samovar_ungz_")
    os.close(fd)
    tmp_path = Path(tmp)
    try:
        gunzip_file(src, tmp_path, remove_source=False)
        yield str(tmp_path)
    finally:
        try:
            tmp_path.unlink()
        except OSError:
            pass


def sibling_candidates(path: PathLike, taxid: str, suffixes: Sequence[str]) -> List[Path]:
    """Candidate sibling files for a genome (plain and gzipped suffixes)."""
    p = as_path(path)
    parent = p.parent
    stem = sequence_stem(p)
    names = []
    for suffix in suffixes:
        names.extend(
            [
                parent / f"{taxid}{suffix}",
                parent / f"{taxid}{suffix}.gz" if not suffix.endswith(".gz") else parent / f"{taxid}{suffix}",
                parent / f"{stem}{suffix}",
                parent / f"{stem}{suffix}.gz" if not suffix.endswith(".gz") else parent / f"{stem}{suffix}",
            ]
        )
    # Also try replacing the current extension.
    uncompressed = Path(strip_compression_suffix(str(p)))
    for suffix in suffixes:
        names.append(uncompressed.with_suffix(suffix if suffix.startswith(".") else f".{suffix}"))
        if not suffix.endswith(".gz"):
            names.append(Path(str(uncompressed.with_suffix(suffix)) + ".gz"))
    seen = set()
    ordered: List[Path] = []
    for cand in names:
        key = str(cand)
        if key in seen:
            continue
        seen.add(key)
        ordered.append(cand)
    return ordered
