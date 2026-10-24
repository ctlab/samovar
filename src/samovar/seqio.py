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

# ISS writes ``prefix_R1.fastq``; Illumina/CASAVA often uses ``prefix_1.fastq``.
# Longer ``_R1`` / ``_R2`` suffixes MUST come first: ``foo_R1.fastq`` also ends
# with ``_1.fastq``.
FASTQ_R1_SUFFIXES = (
    "_R1.fastq.gz",
    "_R1.fq.gz",
    "_R1.fastq",
    "_R1.fq",
    "_1.fastq.gz",
    "_1.fq.gz",
    "_1.fastq",
    "_1.fq",
)
FASTQ_R2_SUFFIXES = (
    "_R2.fastq.gz",
    "_R2.fq.gz",
    "_R2.fastq",
    "_R2.fq",
    "_2.fastq.gz",
    "_2.fq.gz",
    "_2.fastq",
    "_2.fq",
)

# Bare FASTQ names are single-end; explicit mate suffixes take precedence.
FASTQ_SE_SUFFIXES = (".fastq.gz", ".fq.gz", ".fastq", ".fq")

PROCESSED_FASTA_EXTS = (
    ".fa.gz",
    ".fa",
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


def is_gzip_file(path: PathLike) -> bool:
    """True when the path looks gzipped or the file starts with gzip magic."""
    target = as_path(path)
    if is_gzip_path(target):
        return True
    try:
        with open(target, "rb") as raw:
            return raw.read(2) == b"\x1f\x8b"
    except OSError:
        return False


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
    """Numeric taxid or assembly accession from a processed/raw FASTA name."""
    name = strip_compression_suffix(as_path(path).name)
    lower = name.lower()
    for marker in ("-processed.fasta", "-processed.fna", "-processed.fa"):
        if lower.endswith(marker):
            stem = name[: -len(marker)]
            return stem or None
    stem = sequence_stem(path)
    return stem or None


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
    gzipped = is_gzip_path(path)
    if not gzipped and "r" in mode:
        gzipped = is_gzip_file(path)
    if gzipped:
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
    if not is_gzip_file(src):
        return src
    if dest is not None:
        out = as_path(dest)
    elif is_gzip_path(src):
        out = Path(strip_compression_suffix(str(src)))
    else:
        out = src
    out.parent.mkdir(parents=True, exist_ok=True)
    same = src.resolve() == out.resolve()
    write_to = Path(str(out) + ".ungz.tmp") if same else out
    with gzip.open(src, "rb") as fin, open(write_to, "wb") as fout:
        shutil.copyfileobj(fin, fout)
    if same:
        write_to.replace(out)
    elif remove_source:
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
    ident = str(taxid).strip()
    suffix = ".fa.gz" if gzip_genomes else ".fa"
    return as_path(genome_dir) / f"{ident}{suffix}"


def find_existing_processed_genome(genome_dir: PathLike, taxid: str) -> Optional[Path]:
    ident = str(taxid).strip()
    directory = as_path(genome_dir)
    names = [ident]
    search_dirs = [directory, directory / "processed", directory / ".genomes" / "processed"]
    preferred = (
        ".fa.gz",
        "-processed.fasta.gz",
        "-processed.fna.gz",
        "-processed.fa.gz",
        "-processed.fasta",
        "-processed.fna",
        "-processed.fa",
        ".fa",
    )
    for folder in search_dirs:
        if not folder.is_dir() and folder != directory:
            continue
        for name in names:
            for ext in preferred:
                candidate = folder / f"{name}{ext}"
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
    ident = str(taxid).strip()
    found = find_existing_processed_genome(genome_dir, ident)
    if found is not None:
        return str(found)
    directory = as_path(genome_dir)
    for ext in genome_lookup_extensions():
        candidate = directory / f"{ident}{ext}"
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
        return FASTQ_R1_SUFFIXES + FASTQ_SE_SUFFIXES
    if mate == "R2":
        return FASTQ_R2_SUFFIXES
    raise ValueError(f"mate must be R1 or R2, not {mate!r}")


def sample_name_from_r1(path: PathLike) -> str:
    name = as_path(path).name
    if name.endswith(FASTQ_R2_SUFFIXES):
        raise ValueError(f"Not a recognised R1 FASTQ name: {path}")
    for suffix in mate_suffixes("R1"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    raise ValueError(f"Not a recognised R1 FASTQ name: {path}")


def r2_for_r1(r1: PathLike) -> Path:
    path = as_path(r1)
    name = path.name
    for s1, s2 in zip(FASTQ_R1_SUFFIXES, FASTQ_R2_SUFFIXES):
        if name.endswith(s1):
            return path.with_name(name[: -len(s1)] + s2)
    if not name.endswith(FASTQ_R2_SUFFIXES):
        for suffix in FASTQ_SE_SUFFIXES:
            if name.endswith(suffix):
                return path.with_name(name[: -len(suffix)] + "_R2" + suffix)
    raise ValueError(f"Not a recognised R1 FASTQ name: {r1}")


def find_fastq_mate(directory: PathLike, sample: str, mate: str = "R1") -> Optional[Path]:
    directory = as_path(directory)
    candidates = [directory / f"{sample}{suffix}" for suffix in mate_suffixes(mate)]
    found = [path for path in candidates if path.is_file()]
    if len(found) > 1:
        raise ValueError(
            f"Multiple {mate} FASTQ files for sample {sample!r}: "
            + ", ".join(str(path) for path in found)
            + ". Use unique sample names."
        )
    return found[0] if found else None


# CAMISIM harvest writes ``{n}_{tech}_R1.fastq`` then copies to ``{n}_full_``.
# Annotators must only see the ISS-style ``*_full`` names.
_CAMISIM_TECH_TAIL = re.compile(
    r"_(illumina|bgi|mgi|ont|wgsim|nanosim3|art|sequence_type|read_type)$",
    re.IGNORECASE,
)
_CAMISIM_TECH_SAMPLE = re.compile(
    r"^(?P<n>\d+)_(illumina|bgi|mgi|ont|wgsim|nanosim3|art|sequence_type|read_type)$",
    re.IGNORECASE,
)
_FULL_SAMPLE = re.compile(r"^(?P<n>\d+)_full$", re.IGNORECASE)


def _drop_duplicate_camisim_tech_samples(samples: Sequence[str]) -> List[str]:
    stems = set(samples)
    full_ns = {m.group("n") for s in samples if (m := _FULL_SAMPLE.match(s))}
    out: List[str] = []
    for s in samples:
        m = _CAMISIM_TECH_SAMPLE.match(s)
        if m and m.group("n") in full_ns:
            continue
        tail = _CAMISIM_TECH_TAIL.search(s)
        if tail:
            parent = s[: tail.start()]
            if parent in stems or f"{parent}_full" in stems:
                continue
        out.append(s)
    return out


def list_fastq_samples(directory: PathLike, *, include_libraries: bool = False) -> List[str]:
    directory = as_path(directory)
    if not directory.is_dir():
        return []
    samples: List[str] = []
    seen = set()
    for path in sorted(directory.iterdir()):
        if not path.is_file():
            continue
        name = path.name
        if name.endswith(FASTQ_R2_SUFFIXES):
            continue
        for suffix in mate_suffixes("R1"):
            if name.endswith(suffix):
                sample = name[: -len(suffix)]
                if sample not in seen:
                    seen.add(sample)
                    samples.append(sample)
                break
    for sample in samples:
        for mate in ("R1", "R2"):
            find_fastq_mate(directory, sample, mate)
    return samples if include_libraries else _drop_duplicate_camisim_tech_samples(samples)


def list_r1_files(directory: PathLike, *, include_libraries: bool = False) -> List[Path]:
    directory = as_path(directory)
    files: List[Path] = []
    for sample in list_fastq_samples(directory, include_libraries=include_libraries):
        found = find_fastq_mate(directory, sample, "R1")
        if found is not None:
            files.append(found)
    return files


def fastq_pair_paths(prefix: PathLike, gzip_reads: bool = False) -> Tuple[str, str]:
    """Return ``{prefix}_R1.fastq[.gz]`` and R2 paths (not necessarily existing)."""
    prefix_s = str(prefix)
    ext = ".fastq.gz" if gzip_reads else ".fastq"
    return f"{prefix_s}_R1{ext}", f"{prefix_s}_R2{ext}"


def fastq_record_id(header: str) -> str:
    """Read id with only a trailing ``/1`` or ``/2`` mate suffix removed.

    ``r/10`` stays ``r/10``. A suffix in the middle of the token is not a mate mark.
    """
    token = str(header or "").strip()
    if token.startswith("@"):
        token = token[1:]
    token = token.split()[0] if token else ""
    if not token:
        raise ValueError(f"FASTQ record has an empty read id: {header!r}")
    if len(token) >= 2 and token[-2] == "/" and token[-1] in "12":
        return token[:-2]
    return token


def iter_fastq_records(path: PathLike) -> Iterator[Tuple[str, str, str, str]]:
    """Yield ``(header, seq, plus, qual)`` lines, including their newlines.

    A blank line or a record cut off before the quality line is an error.
    """
    target = as_path(path)
    with open_text(target) as handle:
        line_no = 0
        while True:
            header = handle.readline()
            if not header:
                return
            line_no += 1
            if not header.strip():
                raise ValueError(f"Blank line in FASTQ {target} at line {line_no}")
            if not header.startswith("@"):
                raise ValueError(f"FASTQ record at {target}:{line_no} does not start with @")
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            line_no += 3
            if not (seq and plus and qual):
                raise ValueError(f"Truncated FASTQ record in {target} at line {line_no}")
            if not plus.startswith("+"):
                raise ValueError(f"FASTQ record at {target}:{line_no} is missing the + line")
            yield header, seq, plus, qual


DEFAULT_MIN_FASTQ_BYTES = 1
DEFAULT_MIN_FASTQ_READS = 1


def fastq_has_reads(
    path: PathLike,
    *,
    min_bytes: Optional[int] = None,
    min_reads: Optional[int] = None,
) -> bool:
    """True when ``path`` is a FASTQ with enough bytes and complete records.

    Empty or 1-byte stubs (and gzip files that decompress to no records) fail.
    Tiny valid test records such as ``@a\\nA\\n+\\nI\\n`` still pass.
    """
    target = as_path(path)
    try:
        if not target.is_file():
            return False
        size = target.stat().st_size
    except OSError:
        return False
    if size <= 0:
        return False
    need_bytes = DEFAULT_MIN_FASTQ_BYTES if min_bytes is None else int(min_bytes)
    if need_bytes > 0 and size < need_bytes:
        return False
    need_reads = DEFAULT_MIN_FASTQ_READS if min_reads is None else int(min_reads)
    if need_reads <= 0:
        return True
    try:
        n = 0
        for _ in iter_fastq_records(target):
            n += 1
            if n >= need_reads:
                return True
    except (OSError, ValueError):
        return False
    return False


def require_generated_samples(directory: PathLike, samples: Sequence[str], annotator: str = "full") -> None:
    """Reject a successful simulation that left any requested sample empty."""
    if not samples:
        raise ValueError("No samples requested for read generation")
    empty = []
    for sample in samples:
        path = find_fastq_mate(directory, f"{sample}_{annotator}", "R1")
        if path is None or not fastq_has_reads(path):
            empty.append(str(sample))
    if empty:
        raise RuntimeError(
            f"Read generation produced empty sample(s) {', '.join(empty)} under {directory}. "
            "Increase --total_reads or check the simulator output."
        )


def has_r1_reads(
    directory: PathLike,
    *,
    min_bytes: Optional[int] = None,
    min_reads: Optional[int] = None,
) -> bool:
    """True when the directory has at least one R1 FASTQ with real records."""
    for r1 in list_r1_files(directory):
        if fastq_has_reads(r1, min_bytes=min_bytes, min_reads=min_reads):
            return True
    return False


def write_text_lines(path: PathLike, chunks: Iterable[str]) -> None:
    dest = as_path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with open_text(dest, "wt") as handle:
        for chunk in chunks:
            handle.write(chunk)


def concat_fastq_files(sources: Sequence[PathLike], dest: PathLike) -> Path:
    dest_p = as_path(dest)
    dest_p.parent.mkdir(parents=True, exist_ok=True)
    with open_text(dest_p, "wt") as out:
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
    return dest_p


def concat_r1_fastqs(directory: PathLike, dest: PathLike) -> Path:
    return concat_fastq_files(list_r1_files(directory), dest)


def _fastq_r2_name(name: str) -> bool:
    return any(name.endswith(suffix) for suffix in FASTQ_R2_SUFFIXES)


def _stem_before_r2(name: str) -> Optional[str]:
    for suffix in FASTQ_R2_SUFFIXES:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return None


def single_end_sample_names(directory: PathLike) -> List[str]:
    """Sample stems in ``directory`` that have R1 and no non-empty R2."""
    directory = as_path(directory)
    names: List[str] = []
    for r1 in list_r1_files(directory):
        sample = sample_name_from_r1(r1)
        mate = find_fastq_mate(directory, sample, "R2")
        if mate is not None and mate.is_file() and mate.stat().st_size > 0:
            continue
        if sample not in names:
            names.append(sample)
    return names


def canonicalize_single_end_fastqs(directory: PathLike) -> dict[Path, Path]:
    """Use bare FASTQ names for generated unpaired libraries; preserve pairs."""
    renamed = {}
    for r1 in list_r1_files(directory, include_libraries=True):
        sample = sample_name_from_r1(r1)
        r2 = find_fastq_mate(directory, sample, "R2")
        if r2 is not None and r2.stat().st_size > 0:
            continue
        for suffix in FASTQ_R1_SUFFIXES:
            if r1.name.endswith(suffix):
                dest = r1.with_name(sample + suffix[suffix.index("."):])
                r1.rename(dest)
                renamed[r1] = dest
                break
        if r2 is not None:
            r2.unlink()
    return renamed


def drop_r2_if_source_single_end(output_dir: PathLike) -> List[Path]:
    """Delete regenerated R2 when the matching input sample was single-end.

    ISS and ONT simulators write ``*_R2.fastq`` even for an unpaired source.
    A present R2 file, including an empty one, makes later stages call the
    paired CLI. The surviving unpaired file uses a bare FASTQ name; its record
    count is the pair count the abundance table already uses.
    """
    dest = as_path(output_dir)
    if not dest.is_dir():
        return []
    root = dest.resolve().parent
    source: Optional[Path] = None
    for name in ("initial_trimmed", "initial"):
        candidate = root / name
        if candidate.is_dir() and candidate.resolve() != dest.resolve() and list_r1_files(candidate):
            source = candidate
            break
    se_names = set(single_end_sample_names(source)) if source is not None else set()
    source_names = list_fastq_samples(source) if source is not None else []
    all_se = bool(se_names) and len(se_names) == len(source_names)
    removed: List[Path] = []
    for path in list(dest.iterdir()):
        if not path.is_file() or not _fastq_r2_name(path.name):
            continue
        stem = _stem_before_r2(path.name) or ""
        # A hybrid canonical R1 is SE, but its Illumina library is still paired.
        tail = _CAMISIM_TECH_TAIL.search(stem)
        if tail and tail.group(1).lower() in {"illumina", "art", "wgsim", "bgi", "mgi"}:
            base = stem[:tail.start()]
            if any(find_fastq_mate(dest, f"{base}_{tech}") for tech in ("ont", "nanosim3")):
                continue
        if all_se:
            drop = True
        else:
            hits = [
                sample
                for sample in source_names
                if stem == sample or stem.startswith(sample + "_")
            ]
            owner = max(hits, key=len) if hits else None
            drop = owner in se_names
        if not drop:
            continue
        try:
            path.unlink()
        except OSError:
            continue
        removed.append(path)
    canonicalize_single_end_fastqs(dest)
    return removed


def link_or_copy_reads(src_dir: PathLike, dst_dir: PathLike, *, validate_only: bool = False) -> int:
    """Validate sample identities, then symlink/copy pairs unless validate_only."""
    src = as_path(src_dir).resolve()
    dst = as_path(dst_dir).resolve()
    sources = list_r1_files(src, include_libraries=True)
    list_fastq_samples(dst)
    if src == dst:
        return 0
    import filecmp

    filecmp.clear_cache()  # Copies must be rechecked even if file timestamps did not change.
    transfers = []
    for r1 in sources:
        sample = sample_name_from_r1(r1)
        for mate in ("R1", "R2"):
            src_file = find_fastq_mate(src, sample, mate)
            existing = find_fastq_mate(dst, sample, mate)
            if src_file is not None and existing is not None:
                same = existing.name == src_file.name and (
                    existing.samefile(src_file) or (
                        not existing.is_symlink() and filecmp.cmp(existing, src_file, shallow=False)
                    )
                )
                if not same:
                    raise ValueError(f"Sample name collision for {sample!r} ({mate}): {src_file} and {existing}")
            if src_file is not None:
                transfers.append(src_file)
    if validate_only:
        return 0
    dst.mkdir(parents=True, exist_ok=True)
    n = 0
    for src_file in transfers:
        target = dst / src_file.name
        if target.exists() or target.is_symlink():
            continue
        try:
            target.symlink_to(src_file.resolve())
        except OSError:
            shutil.copy2(src_file, target)
        n += 1
    return n


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


@contextmanager
def assembly_read_libraries(r1: PathLike, r2: PathLike = ""):
    """Yield (short R1, short R2, ONT, typed) without guessing lost pairing.

    Generated hybrids retain libraries beside the canonical classifier input.
    Older flattened inputs can only recover short reads as an unpaired library.
    Untagged inputs retain their caller's existing interpretation.
    """
    r1 = as_path(r1)
    try:
        stem = sample_name_from_r1(r1)
    except ValueError:
        stem = r1.stem
    bases = [stem]
    if stem.endswith("_full"):
        bases.append(stem[:-5])
    for base in bases:
        short = []
        long = []
        present = set()
        for tech in ("illumina", "art", "wgsim", "bgi", "mgi", "ont", "nanosim3"):
            path = find_fastq_mate(r1.parent, f"{base}_{tech}")
            kind = "long" if tech in {"ont", "nanosim3"} else "short"
            if path:
                present.add(kind)
                if path.stat().st_size:
                    (long if kind == "long" else short).append(path)
        # A lone alias is not enough to replace a potentially mixed canonical R1.
        if present != {"short", "long"}:
            continue
        if len(short) > 1 or len(long) > 1:
            raise ValueError("Hybrid assembly requires one short-read and one ONT library")
        mate = find_fastq_mate(r1.parent, sample_name_from_r1(short[0]), "R2") if short else None
        yield str(short[0]) if short else "", str(mate or ""), str(long[0]) if long else "", True
        return
    with tempfile.TemporaryDirectory(prefix="samovar_libraries_") as tmp:
        short_path = Path(tmp) / "short.fastq"
        long_path = Path(tmp) / "ont.fastq"
        tagged = untagged = 0
        with short_path.open("w") as sh, long_path.open("w") as lg:
            for path in (r1, as_path(r2) if r2 else None):
                if path is None or not path.is_file():
                    continue
                for record in iter_fastq_records(path):
                    match = re.search(r"(?:^|[|\s])read_type:([^|\s/]+)", record[0].lstrip("@"))
                    if not match:
                        untagged += 1
                        continue
                    tech = match.group(1).lower()
                    if tech not in {"illumina", "art", "wgsim", "bgi", "mgi", "ont", "nanosim", "nanosim3"}:
                        raise ValueError(f"Unsupported assembly read_type: {tech}")
                    tagged += 1
                    (lg if tech in {"ont", "nanosim", "nanosim3"} else sh).writelines(record)
        if tagged and untagged:
            raise ValueError("Cannot assemble a mixture of tagged and untagged read types")
        if not tagged:
            yield str(r1), str(r2 or ""), "", False
        elif long_path.stat().st_size:
            yield str(short_path) if short_path.stat().st_size else "", "", str(long_path), True
        else:
            # Pure short reads retain their original pairing.
            yield str(r1), str(r2 or ""), "", True
