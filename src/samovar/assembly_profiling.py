"""Assembly-based annotation: nested contracts plus a composite annotator.

Outer CLI matches ``annotator`` (``-i/-I/-d/-o/-t`` → seq/taxID TSV). Inner
stages are importable groups (assembler, gene_caller, binner, binner_qc,
binner_combine, mag_taxonomy, aligner, mag_quantifier, taxon_quantifier,
read_assigner).
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

from samovar.paths import resolve_executable

MAG_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz")
BINARY_NAMES = {
    "megahit": "megahit",
    "minimap2": "minimap2",
    "samtools": "samtools",
    "coverm": "coverm",
    "checkm2": "checkm2",
    "gtdbtk": "gtdbtk",
    "prodigal": "prodigal",
    "dastool": "DAS_Tool",
    "DAS_Tool": "DAS_Tool",
    "anvio": "anvi-gen-contigs-database",
    "anvi-gen-contigs-database": "anvi-gen-contigs-database",
    "anvi-profile": "anvi-profile",
    "anvi-merge": "anvi-merge",
    "anvi-cluster-contigs": "anvi-cluster-contigs",
    "anvi-summarize": "anvi-summarize",
    "metabat2": "metabat2",
}
STAGE_NAMES = (
    "assembler",
    "gene_caller",
    "binner",
    "binner_qc",
    "binner_combine",
    "mag_taxonomy",
    "aligner",
    "mag_quantifier",
    "taxon_quantifier",
    "read_assigner",
    "annotator",
)


def as_path(value: Any) -> Path:
    return Path(str(value)).expanduser()


def which_tool(name: str, default: str = "") -> str:
    binary = BINARY_NAMES.get(name, name)
    raw = default or binary
    resolved = resolve_executable(raw, tool_key=name) or resolve_executable(
        binary, tool_key=binary
    )
    first = str(resolved).split()[0] if resolved else ""
    if first and Path(first).expanduser().is_file() and os.access(first, os.X_OK):
        return first
    from samovar.paths import sidecar_envs_dir

    envs = sidecar_envs_dir()
    candidates = []
    if envs.is_dir():
        candidates.extend(envs.glob(f"*/bin/{binary}"))
    for folder in (
        name,
        Path(binary).name.split("-")[0],
        "anvio",
        "dastool",
        "gtdbtk",
        "checkm2",
        "minimap2",
        "megahit",
        "prodigal",
    ):
        candidates.append(envs / folder / "bin" / binary)
    for cand in candidates:
        if cand.is_file() and os.access(cand, os.X_OK):
            return str(cand)
    found = shutil.which(binary) or shutil.which(name) or shutil.which(raw)
    if found:
        return found
    raise FileNotFoundError(
        f"{name} is not installed. Import it with samovar tools import "
        f"or ./install.sh {name}"
    )


def gtdbtk_db_ready(path: str) -> bool:
    """True when GTDB-Tk data is extracted (not only a download directory)."""
    root = as_path(path)
    if not root.is_dir():
        return False
    if (root / "metadata").is_dir() or (root / "taxonomy").is_dir():
        return True
    if (root / "tigrfam").is_dir() or (root / "hmm_metadata").is_dir():
        return True
    return False


def gtdbtk_data_path(explicit: str = "") -> str:
    if explicit and gtdbtk_db_ready(explicit):
        return str(as_path(explicit))
    env = os.environ.get("GTDBTK_DATA_PATH", "").strip()
    if env and gtdbtk_db_ready(env):
        return env
    for shared in (
        Path("/mnt/tank/scratch/partition-metagenomics/databases/GTDB"),
        Path("/mnt/tank/scratch/partition-metagenomics/databases/gtdbtk"),
    ):
        if gtdbtk_db_ready(str(shared)):
            return str(shared)
    return explicit or env


def assembly_slot_extra(
    *,
    assembler: str = "",
    gene_caller: str = "",
    binners: Optional[Sequence[str]] = None,
    binner_qc: str = "",
    binner_combine: str = "",
    mag_taxonomy: str = "",
    aligner: str = "",
    mag_quantifier: str = "",
) -> str:
    parts: List[str] = []
    if assembler:
        parts.extend(["--assembler", str(assembler)])
    if gene_caller:
        parts.extend(["--gene-caller", str(gene_caller)])
    for name in binners or []:
        if name:
            parts.extend(["--binner", str(name)])
    if binner_qc:
        parts.extend(["--binner-qc", str(binner_qc)])
    if binner_combine:
        parts.extend(["--binner-combine", str(binner_combine)])
    if mag_taxonomy:
        parts.extend(["--mag-taxonomy", str(mag_taxonomy)])
    if aligner:
        parts.extend(["--aligner", str(aligner)])
    if mag_quantifier:
        parts.extend(["--mag-quantifier", str(mag_quantifier)])
    return " ".join(parts)


def parse_slot_extra(extra: str) -> Dict[str, Any]:
    tokens = str(extra or "").split()
    out: Dict[str, Any] = {"binners": []}
    i = 0
    while i < len(tokens):
        tok = tokens[i]
        nxt = tokens[i + 1] if i + 1 < len(tokens) else ""
        if tok == "--assembler" and nxt:
            out["assembler"] = nxt
            i += 2
            continue
        if tok in {"--gene-caller", "--gene_caller", "--gene-annotation"} and nxt:
            out["gene_caller"] = nxt
            i += 2
            continue
        if tok == "--binner" and nxt:
            out["binners"].append(nxt)
            i += 2
            continue
        if tok == "--binner-qc" and nxt:
            out["binner_qc"] = nxt
            i += 2
            continue
        if tok == "--binner-combine" and nxt:
            out["binner_combine"] = nxt
            i += 2
            continue
        if tok == "--mag-taxonomy" and nxt:
            out["mag_taxonomy"] = nxt
            i += 2
            continue
        if tok == "--aligner" and nxt:
            out["aligner"] = nxt
            i += 2
            continue
        if tok == "--mag-quantifier" and nxt:
            out["mag_quantifier"] = nxt
            i += 2
            continue
        i += 1
    return out


def _run(cmd: Sequence[str], **kwargs: Any) -> None:
    env = dict(kwargs.pop("env", None) or os.environ)
    env.setdefault("PYTHONNOUSERSITE", "1")
    subprocess.check_call([str(c) for c in cmd], env=env, **kwargs)


def _is_mag_fasta(path: Path) -> bool:
    name = path.name.lower()
    return any(name.endswith(suf) for suf in MAG_SUFFIXES)


def _fasta_has_records(path: Path) -> bool:
    src = as_path(path)
    if not src.is_file() or src.stat().st_size == 0:
        return False
    with src.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                return True
    return False


def _iter_mag_fastas(mag_dir: Path) -> List[Path]:
    if mag_dir.is_file() and _is_mag_fasta(mag_dir):
        return [mag_dir]
    if not mag_dir.is_dir():
        return []
    found: List[Path] = []
    for path in sorted(mag_dir.rglob("*")):
        if not path.is_file() or not _is_mag_fasta(path):
            continue
        rel = path.relative_to(mag_dir)
        if any(part.startswith(".") for part in rel.parts):
            continue
        found.append(path)
    mag_bins = [p for p in found if p.stem.lower().startswith(("mag", "bin"))]
    if mag_bins and any(p.stem.lower() == "contigs" for p in found):
        return mag_bins
    return found


def _prefix_fasta_headers(src: Path, dest: Path, mag_id: str) -> Path:
    """Rewrite ``>c1`` to ``>mag_id~c1`` so CoverM/minimap contig names are unique."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    parts: List[str] = []
    prefix = f"{mag_id}~"
    try:
        text = src.read_text(encoding="utf-8", errors="replace")
    except OSError:
        dest.write_text("", encoding="utf-8")
        return dest
    for line in text.splitlines():
        if line.startswith(">"):
            contig = line[1:].split()[0]
            if contig.startswith(prefix):
                parts.append(f">{contig}")
            else:
                parts.append(f">{prefix}{contig}")
        else:
            parts.append(line)
    dest.write_text("\n".join(parts) + ("\n" if parts else ""), encoding="utf-8")
    return dest


def _copy_fasta(src: Path, dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)


def _touch_empty_fasta(dest: Path) -> None:
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("", encoding="utf-8")


def run_megahit(
    r1: str,
    r2: str = "",
    dest: str = "",
    threads: int = 1,
    extra: Optional[Sequence[str]] = None,
    exe: str = "",
) -> Path:
    """Assemble paired (or single) FASTQ with MegaHIT; write contig FASTA to dest."""
    out = as_path(dest)
    work = out.parent / (out.stem + ".megahit")
    if work.exists():
        shutil.rmtree(work)
    cmd = [which_tool("megahit", exe), "-t", str(int(threads)), "-o", str(work), "-f"]
    if r2 and as_path(r2).is_file() and as_path(r2).stat().st_size > 0:
        cmd.extend(["-1", str(r1), "-2", str(r2)])
    else:
        cmd.extend(["-r", str(r1)])
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    contigs = work / "final.contigs.fa"
    out.parent.mkdir(parents=True, exist_ok=True)
    if contigs.is_file() and contigs.stat().st_size > 0:
        shutil.copy2(contigs, out)
    else:
        _touch_empty_fasta(out)
    return out


def run_assembler(
    r1: str,
    r2: str = "",
    dest: str = "",
    threads: int = 1,
    name: str = "megahit",
    extra: Optional[Sequence[str]] = None,
) -> Path:
    key = str(name or "megahit").lower()
    if key in {"megahit", "assembler", ""}:
        return run_megahit(r1, r2, dest, threads=threads, extra=extra)
    exe = which_tool(name)
    out = as_path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "-i", str(r1), "-I", str(r2 or ""), "-o", str(out), "-t", str(int(threads))]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    return out


def _fasta_inputs(src: Path) -> List[Path]:
    if src.is_file() and _is_mag_fasta(src):
        return [src]
    return _iter_mag_fastas(src)


def run_prodigal(
    contigs: str,
    dest: str,
    threads: int = 1,
    extra: Optional[Sequence[str]] = None,
    exe: str = "",
) -> Path:
    """Call genes with Prodigal; write ``*.faa`` and ``*.gff`` under dest."""
    del threads  # Prodigal is single-threaded
    out = as_path(dest)
    out.mkdir(parents=True, exist_ok=True)
    prodigal = which_tool("prodigal", exe)
    files = _fasta_inputs(as_path(contigs))
    if not files and as_path(contigs).is_file():
        files = [as_path(contigs)]
    for fasta in files:
        stem = fasta.stem
        faa = out / f"{stem}.faa"
        gff = out / f"{stem}.gff"
        cmd = [
            prodigal,
            "-i",
            str(fasta),
            "-a",
            str(faa),
            "-o",
            str(gff),
            "-f",
            "gff",
            "-p",
            "meta",
        ]
        cmd.extend(str(x) for x in (extra or []))
        _run(cmd)
    return out


def run_gene_caller(
    contigs: str,
    dest: str,
    threads: int = 1,
    name: str = "prodigal",
    extra: Optional[Sequence[str]] = None,
) -> Path:
    key = str(name or "prodigal").lower().replace("-", "_")
    cand = as_path(name) if name else None
    if cand and cand.is_file():
        out = as_path(dest)
        out.mkdir(parents=True, exist_ok=True)
        exe = str(cand)
        cmd = [exe, "-c", str(contigs), "-o", str(out), "-t", str(int(threads))]
        if exe.endswith(".py"):
            cmd = [sys.executable, *cmd]
        cmd.extend(str(x) for x in (extra or []))
        _run(cmd)
        return out
    if key in {"prodigal", "gene_caller", "gene_annotation", ""}:
        return run_prodigal(contigs, dest, threads=threads, extra=extra)
    exe = which_tool(name)
    out = as_path(dest)
    out.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "-c", str(contigs), "-o", str(out), "-t", str(int(threads))]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    return out


def _write_contig2bin(mag_dir: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for fasta in _iter_mag_fastas(mag_dir):
        mag_id = fasta.stem
        header = ""
        with fasta.open() as handle:
            for line in handle:
                if line.startswith(">"):
                    header = line[1:].split()[0].strip()
                    if header:
                        rows.append((header, mag_id))
    dest.write_text("".join(f"{c}\t{b}\n" for c, b in rows), encoding="utf-8")
    return dest


def _write_identity_mag(contigs_fa: Path, mag_dir: Path) -> Path:
    """One MAG = the assembly (used when clustering is skipped)."""
    mag_dir.mkdir(parents=True, exist_ok=True)
    _copy_fasta(as_path(contigs_fa), mag_dir / "mag1.fa")
    return mag_dir


def run_anvio_binner(
    contigs: str,
    dest: str,
    r1: str = "",
    r2: str = "",
    threads: int = 1,
    extra: Optional[Sequence[str]] = None,
    additional_reads: Optional[Sequence[Tuple[str, str]]] = None,
    skip_clustering: Optional[bool] = None,
) -> Path:
    """Bin contigs with anvi'o; clustering is optional for a single sample.

    Internal contig mapping still runs. ``anvi-cluster-contigs`` needs a merged
    profile (two or more samples). With one sample, clustering is skipped and the
    assembly is written as ``mag1.fa``. Extra FASTQ pairs in ``additional_reads``
    are profiled, merged, and clustered.
    """
    mag_dir = as_path(dest)
    mag_dir.mkdir(parents=True, exist_ok=True)
    work = mag_dir / ".anvio"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    simple = work / "contigs.fa"
    i = 0
    lines = []
    for line in as_path(contigs).read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            i += 1
            lines.append(f">c{i}")
        else:
            lines.append(line)
    simple.write_text("\n".join(lines) + "\n", encoding="utf-8")
    if not _fasta_has_records(simple):
        return _write_identity_mag(simple, mag_dir)
    db = work / "CONTIGS.db"
    gen = which_tool("anvi-gen-contigs-database")
    gen_cmd = [
        gen,
        "-f",
        str(simple),
        "-o",
        str(db),
        "-n",
        "samovar",
        "-T",
        str(int(threads)),
    ]
    env = os.environ.copy()
    try:
        prodigal = which_tool("prodigal")
        env["PATH"] = str(Path(prodigal).parent) + os.pathsep + env.get("PATH", "")
        _run(gen_cmd, env=env)
    except (FileNotFoundError, subprocess.CalledProcessError):
        if db.exists():
            db.unlink()
        try:
            _run(gen_cmd + ["--skip-gene-calling"])
        except (FileNotFoundError, subprocess.CalledProcessError, OSError):
            return _write_identity_mag(simple, mag_dir)
    bam = work / "contigs.bam"
    pairs: List[Tuple[str, str]] = []
    if r1 and as_path(r1).is_file():
        pairs.append((str(r1), str(r2 or "")))
    for extra_pair in additional_reads or []:
        if extra_pair and extra_pair[0] and as_path(extra_pair[0]).is_file():
            pairs.append((str(extra_pair[0]), str(extra_pair[1] if len(extra_pair) > 1 else "")))
    profile_dbs: List[Path] = []
    for i, (pr1, pr2) in enumerate(pairs, start=1):
        sample_bam = work / f"sample{i}.bam" if len(pairs) > 1 else bam
        try:
            run_minimap2(pr1, pr2, str(simple), str(sample_bam), threads=threads)
        except FileNotFoundError:
            continue
        if not sample_bam.is_file() or sample_bam.stat().st_size <= 0:
            continue
        pdir = work / (f"profile_{i}" if len(pairs) > 1 else "profile")
        try:
            _run(
                [
                    which_tool("anvi-profile"),
                    "-i",
                    str(sample_bam),
                    "-c",
                    str(db),
                    "-o",
                    str(pdir),
                    "-S",
                    f"S{i}",
                    "-T",
                    str(int(threads)),
                    "--skip-hierarchical-clustering",
                ]
            )
        except (FileNotFoundError, subprocess.CalledProcessError):
            continue
        pdb = pdir / "PROFILE.db"
        if pdb.is_file():
            profile_dbs.append(pdb)
    want_cluster = skip_clustering is False or (
        skip_clustering is None
        and not any(str(x) in {"--skip-clustering", "--no-cluster"} for x in (extra or []))
    )
    clustered = False
    profile = profile_dbs[0] if profile_dbs else work / "profile" / "PROFILE.db"
    extra_list = [str(x) for x in (extra or []) if str(x) not in {"--skip-clustering", "--no-cluster"}]
    if want_cluster and len(profile_dbs) >= 2:
        merge_dir = work / "merged"
        try:
            merge_cmd = [which_tool("anvi-merge"), "-i"]
            merge_cmd.extend(str(p) for p in profile_dbs)
            merge_cmd.extend(
                [
                    "-o",
                    str(merge_dir),
                    "-c",
                    str(db),
                    "--skip-hierarchical-clustering",
                ]
            )
            merge_env = os.environ.copy()
            merge_bin = Path(merge_cmd[0]).parent
            merge_env["PATH"] = str(merge_bin) + os.pathsep + merge_env.get("PATH", "")
            _run(merge_cmd, env=merge_env)
            merged = merge_dir / "PROFILE.db"
            if merged.is_file():
                profile = merged
                cluster = which_tool("anvi-cluster-contigs")
                cmd = [
                    cluster,
                    "-c",
                    str(db),
                    "-p",
                    str(profile),
                    "--driver",
                    "concoct",
                    "-C",
                    "BINS",
                    "-T",
                    str(int(threads)),
                    "--just-do-it",
                ]
                cmd.extend(extra_list)
                _run(cmd, env=merge_env)
                clustered = True
        except (FileNotFoundError, subprocess.CalledProcessError):
            clustered = False
    elif want_cluster and len(profile_dbs) == 1:
        clustered = False
    if not clustered:
        return _write_identity_mag(simple, mag_dir)
    summary = work / "summary"
    try:
        _run(
            [
                which_tool("anvi-summarize"),
                "-c",
                str(db),
                "-p",
                str(profile),
                "-C",
                "BINS",
                "-o",
                str(summary),
            ]
        )
    except Exception:
        pass
    bins_fasta = []
    if summary.is_dir():
        bins_fasta = [p for p in summary.rglob("*") if p.is_file() and _is_mag_fasta(p)]
    if not bins_fasta:
        return _write_identity_mag(simple, mag_dir)
    for i, fasta in enumerate(bins_fasta, start=1):
        _copy_fasta(fasta, mag_dir / f"mag{i}.fa")
    return mag_dir


def run_binner(
    contigs: str,
    dest: str,
    r1: str = "",
    r2: str = "",
    threads: int = 1,
    name: str = "anvio",
    extra: Optional[Sequence[str]] = None,
    additional_reads: Optional[Sequence[Tuple[str, str]]] = None,
) -> Path:
    key = str(name or "anvio").lower().replace("-", "_")
    if key in {"anvio", "anvi", "anvio_binner"}:
        return run_anvio_binner(
            contigs,
            dest,
            r1=r1,
            r2=r2,
            threads=threads,
            extra=extra,
            additional_reads=additional_reads,
        )
    if key in {"metabat2", "metabat"}:
        exe = which_tool("metabat2")
        mag_dir = as_path(dest)
        mag_dir.mkdir(parents=True, exist_ok=True)
        prefix = mag_dir / "bin"
        cmd = [exe, "-i", str(contigs), "-o", str(prefix), "-t", str(int(threads))]
        cmd.extend(str(x) for x in (extra or []))
        _run(cmd)
        return mag_dir
    exe = which_tool(name)
    mag_dir = as_path(dest)
    mag_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        exe,
        "-c",
        str(contigs),
        "-i",
        str(r1 or ""),
        "-I",
        str(r2 or ""),
        "-o",
        str(mag_dir),
        "-t",
        str(int(threads)),
    ]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    return mag_dir


def run_checkm2(
    mag_dir: str,
    dest: str,
    db: str = "",
    threads: int = 1,
    extra: Optional[Sequence[str]] = None,
) -> Path:
    out = as_path(dest)
    work = out.parent / (out.stem + ".checkm2")
    work.mkdir(parents=True, exist_ok=True)
    cmd = [
        which_tool("checkm2"),
        "predict",
        "--input",
        str(mag_dir),
        "--output-directory",
        str(work),
        "--threads",
        str(int(threads)),
        "--force",
    ]
    if db:
        cmd.extend(["--database_path", str(db)])
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    report = work / "quality_report.tsv"
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = ["mag_id\tcompleteness\tcontamination\tscore"]
    if report.is_file():
        with report.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for rec in reader:
                name = rec.get("Name") or rec.get("name") or rec.get("Bin Id") or ""
                comp = rec.get("Completeness") or rec.get("completeness") or "0"
                cont = rec.get("Contamination") or rec.get("contamination") or "0"
                try:
                    score = float(comp) - 5.0 * float(cont)
                except ValueError:
                    score = 0.0
                rows.append(f"{name}\t{comp}\t{cont}\t{score}")
    elif _iter_mag_fastas(as_path(mag_dir)):
        for fasta in _iter_mag_fastas(as_path(mag_dir)):
            rows.append(f"{fasta.stem}\t0\t0\t0")
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return out


def run_binner_qc(
    mag_dir: str,
    dest: str,
    db: str = "",
    threads: int = 1,
    name: str = "checkm2",
    extra: Optional[Sequence[str]] = None,
) -> Path:
    key = str(name or "checkm2").lower().replace("-", "_")
    if key in {"checkm2", "binner_qc", ""}:
        return run_checkm2(mag_dir, dest, db=db, threads=threads, extra=extra)
    exe = which_tool(name)
    out = as_path(dest)
    cmd = [exe, "-c", str(mag_dir), "-d", str(db or ""), "-o", str(out), "-t", str(int(threads))]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    return out


def run_dastool_combine(
    mag_dirs: Sequence[str],
    dest: str,
    contigs: str = "",
    threads: int = 1,
    extra: Optional[Sequence[str]] = None,
) -> Path:
    mag_dir = as_path(dest)
    mag_dir.mkdir(parents=True, exist_ok=True)
    dirs = [as_path(d) for d in mag_dirs if as_path(d).exists()]
    if not dirs:
        return mag_dir
    if len(dirs) == 1:
        for fasta in _iter_mag_fastas(dirs[0]):
            _copy_fasta(fasta, mag_dir / fasta.name)
        return mag_dir
    work = mag_dir / ".dastool"
    work.mkdir(parents=True, exist_ok=True)
    tables = []
    labels = []
    for i, folder in enumerate(dirs):
        table = work / f"bins{i}.tsv"
        _write_contig2bin(folder, table)
        tables.append(str(table))
        labels.append(folder.name or f"binner{i}")
    contig_fa = contigs
    if not contig_fa:
        fasta = _iter_mag_fastas(dirs[0])
        contig_fa = str(fasta[0]) if fasta else ""
    cmd = [
        which_tool("DAS_Tool", "dastool"),
        "-i",
        ",".join(tables),
        "-l",
        ",".join(labels),
        "-c",
        str(contig_fa),
        "-o",
        str(work / "das"),
        "--write_bins",
        "-t",
        str(int(threads)),
    ]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    bins_out = work / "das_DASTool_bins"
    if not bins_out.is_dir():
        bins_out = work
    found = list(work.rglob("*.fa")) + list(work.rglob("*.fna"))
    if not found:
        for folder in dirs:
            for fasta in _iter_mag_fastas(folder):
                _copy_fasta(fasta, mag_dir / f"{folder.name}_{fasta.name}")
        return mag_dir
    for fasta in found:
        if ".dastool" in str(fasta) and "bins" not in fasta.name.lower() and "bin" not in str(fasta.parent):
            continue
        _copy_fasta(fasta, mag_dir / fasta.name)
    return mag_dir


def run_binner_combine(
    mag_dirs: Sequence[str],
    dest: str,
    qc_table: str = "",
    contigs: str = "",
    threads: int = 1,
    name: str = "dastool",
    extra: Optional[Sequence[str]] = None,
) -> Path:
    key = str(name or "dastool").lower().replace("-", "_")
    if key in {"dastool", "das_tool", "binner_combine", ""}:
        return run_dastool_combine(mag_dirs, dest, contigs=contigs, threads=threads, extra=extra)
    if key in {"identity", "copy", "one"}:
        mag_dir = as_path(dest)
        mag_dir.mkdir(parents=True, exist_ok=True)
        src = as_path(mag_dirs[0]) if mag_dirs else mag_dir
        for fasta in _iter_mag_fastas(src):
            _copy_fasta(fasta, mag_dir / fasta.name)
        return mag_dir
    exe = which_tool(name)
    out = as_path(dest)
    cmd = [
        exe,
        "-c",
        ",".join(str(d) for d in mag_dirs),
        "-q",
        str(qc_table or ""),
        "-o",
        str(out),
        "-t",
        str(int(threads)),
    ]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    return out


def _gtdb_taxid_from_lineage(lineage: str) -> str:
    from samovar.taxonomy import ncbi_taxid_from_gtdb_lineage

    return ncbi_taxid_from_gtdb_lineage(lineage)


def run_gtdbtk(
    mag_dir: str,
    dest: str,
    db: str = "",
    threads: int = 1,
    extra: Optional[Sequence[str]] = None,
) -> Path:
    out = as_path(dest)
    work = out.parent / (out.stem + ".gtdbtk")
    env = os.environ.copy()
    data = gtdbtk_data_path(db)
    if data:
        env["GTDBTK_DATA_PATH"] = str(data)
    if db and gtdbtk_db_ready(db):
        env["GTDBTK_DATA_PATH"] = str(db)
    if not env.get("GTDBTK_DATA_PATH") or not gtdbtk_db_ready(env["GTDBTK_DATA_PATH"]):
        raise FileNotFoundError(
            "GTDB-Tk data is not ready. Extract the package under "
            "/mnt/tank/scratch/partition-metagenomics/databases/GTDB "
            "or set GTDBTK_DATA_PATH."
        )
    try:
        prodigal = which_tool("prodigal")
        env["PATH"] = str(Path(prodigal).parent) + os.pathsep + env.get("PATH", "")
    except FileNotFoundError:
        pass
    if work.exists():
        shutil.rmtree(work)
    cmd = [
        which_tool("gtdbtk"),
        "classify_wf",
        "--genome_dir",
        str(mag_dir),
        "--out_dir",
        str(work),
        "--cpus",
        str(int(threads)),
        "--force",
        "--extension",
        "fa",
    ]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd, env=env)
    rows = ["mag_id\ttaxid\tlineage"]
    for summary in work.rglob("*summary.tsv"):
        with summary.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for rec in reader:
                mag = rec.get("user_genome") or rec.get("Name") or ""
                lineage = rec.get("classification") or rec.get("fastani_taxonomy") or ""
                taxid = _gtdb_taxid_from_lineage(lineage)
                rows.append(f"{mag}\t{taxid}\t{lineage}")
    if len(rows) == 1:
        for fasta in _iter_mag_fastas(as_path(mag_dir)):
            rows.append(f"{fasta.stem}\t0\t")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return out


def run_mag_taxonomy(
    mag_dir: str,
    dest: str,
    db: str = "",
    threads: int = 1,
    name: str = "gtdbtk",
    extra: Optional[Sequence[str]] = None,
) -> Path:
    key = str(name or "gtdbtk").lower().replace("-", "_")
    if key in {"gtdbtk", "gtdb", "mag_taxonomy", ""}:
        return run_gtdbtk(mag_dir, dest, db=db, threads=threads, extra=extra)
    exe = which_tool(name)
    out = as_path(dest)
    cmd = [exe, "-c", str(mag_dir), "-d", str(db or ""), "-o", str(out), "-t", str(int(threads))]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    return out


def contig_to_mag_map(mag_dir: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for fasta in _iter_mag_fastas(mag_dir):
        mag = fasta.stem
        try:
            text = fasta.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for line in text.splitlines():
            if line.startswith(">"):
                contig = line[1:].split()[0]
                out[contig] = mag
                out[f"{mag}~{contig}"] = mag
    return out


def _concat_mag_fasta(mag_dir: Path, dest: Path) -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    parts = []
    for fasta in _iter_mag_fastas(mag_dir):
        mag = fasta.stem
        try:
            text = fasta.read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        for line in text.splitlines():
            if line.startswith(">"):
                contig = line[1:].split()[0]
                if contig.startswith(f"{mag}~"):
                    parts.append(f">{contig}")
                else:
                    parts.append(f">{mag}~{contig}")
            else:
                parts.append(line)
        if parts and parts[-1] != "":
            parts.append("")
    dest.write_text("\n".join(parts) + ("\n" if parts else ""), encoding="utf-8")
    return dest


def run_minimap2(
    r1: str,
    r2: str,
    reference: str,
    dest: str,
    threads: int = 1,
    preset: str = "sr",
    extra: Optional[Sequence[str]] = None,
) -> Path:
    out = as_path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)
    ref = as_path(reference)
    if ref.is_dir():
        ref = _concat_mag_fasta(ref, out.parent / (out.stem + ".mags.fa"))
    mm = which_tool("minimap2")
    samtools = which_tool("samtools")
    cmd = [mm, "-ax", preset, "-t", str(int(threads)), str(ref), str(r1)]
    if r2 and as_path(r2).is_file() and as_path(r2).stat().st_size > 0:
        cmd.append(str(r2))
    cmd.extend(str(x) for x in (extra or []))
    gen = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    sort = subprocess.Popen(
        [samtools, "sort", "-@", str(int(threads)), "-o", str(out)],
        stdin=gen.stdout,
    )
    if gen.stdout:
        gen.stdout.close()
    if sort.wait() != 0 or gen.wait() != 0:
        raise subprocess.CalledProcessError(sort.returncode or gen.returncode, cmd)
    try:
        _run([samtools, "index", str(out)])
    except Exception:
        pass
    return out


def run_aligner(
    r1: str,
    r2: str,
    reference: str,
    dest: str,
    threads: int = 1,
    name: str = "minimap2",
    extra: Optional[Sequence[str]] = None,
) -> Path:
    key = str(name or "minimap2").lower()
    if key in {"minimap2", "aligner", ""}:
        return run_minimap2(r1, r2, reference, dest, threads=threads, extra=extra)
    exe = which_tool(name)
    out = as_path(dest)
    cmd = [
        exe,
        "-i",
        str(r1),
        "-I",
        str(r2 or ""),
        "-r",
        str(reference),
        "-o",
        str(out),
        "-t",
        str(int(threads)),
    ]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    return out


def run_coverm(
    bam: str,
    mag_dir: str,
    dest: str,
    extra: Optional[Sequence[str]] = None,
) -> Path:
    out = as_path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)
    unique_dir = out.parent / (out.stem + ".coverm_genomes")
    if unique_dir.exists():
        shutil.rmtree(unique_dir)
    unique_dir.mkdir(parents=True, exist_ok=True)
    fastas = [
        _prefix_fasta_headers(src, unique_dir / src.name, src.stem)
        for src in _iter_mag_fastas(as_path(mag_dir))
        if _fasta_has_records(src)
    ]
    cmd = [
        which_tool("coverm"),
        "genome",
        "--bam-files",
        str(bam),
        "--methods",
        "count",
        "--min-covered-fraction",
        "0",
    ]
    if fastas:
        cmd.extend(["--genome-fasta-files", *[str(p) for p in fastas]])
    else:
        cmd.extend(
            [
                "--genome-fasta-directory",
                str(unique_dir if unique_dir.is_dir() else mag_dir),
                "--genome-fasta-extension",
                "fa",
            ]
        )
    cmd.extend(str(x) for x in (extra or []))
    text = subprocess.check_output(cmd, text=True)
    lines = [ln for ln in text.splitlines() if ln.strip()]
    rows = ["mag_id\tN"]
    if lines:
        header = lines[0].split("\t")
        mag_i = 0
        n_i = 1 if len(header) > 1 else 1
        for i, name in enumerate(header):
            low = name.lower()
            if "genome" in low or name == "mag_id":
                mag_i = i
            if "count" in low or name == "N":
                n_i = i
        for line in lines[1:]:
            parts = line.split("\t")
            if "unmapped" in parts[0].lower():
                continue
            rows.append(f"{parts[mag_i]}\t{parts[n_i] if n_i < len(parts) else 0}")
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return out


def run_samtools_idxstats(bam: str, mag_dir: str, dest: str) -> Path:
    out = as_path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)
    text = subprocess.check_output([which_tool("samtools"), "idxstats", str(bam)], text=True)
    mapping = contig_to_mag_map(as_path(mag_dir))
    counts: Dict[str, int] = defaultdict(int)
    for line in text.splitlines():
        parts = line.split("\t")
        if len(parts) < 3 or parts[0] == "*":
            continue
        mag = mapping.get(parts[0], parts[0].split("~", 1)[0])
        counts[mag] += int(parts[2])
    if not counts:
        for fasta in _iter_mag_fastas(as_path(mag_dir)):
            counts[fasta.stem] = 0
    rows = ["mag_id\tN"] + [f"{k}\t{v}" for k, v in sorted(counts.items())]
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return out


def run_mag_quantifier(
    bam: str,
    mag_dir: str,
    dest: str,
    name: str = "coverm",
    extra: Optional[Sequence[str]] = None,
) -> Path:
    key = str(name or "coverm").lower()
    if key in {"coverm", "mag_quantifier", ""}:
        try:
            return run_coverm(bam, mag_dir, dest, extra=extra)
        except (FileNotFoundError, subprocess.CalledProcessError):
            return run_samtools_idxstats(bam, mag_dir, dest)
    if key in {"samtools", "idxstats"}:
        return run_samtools_idxstats(bam, mag_dir, dest)
    exe = which_tool(name)
    out = as_path(dest)
    cmd = [exe, "-b", str(bam), "-r", str(mag_dir), "-o", str(out)]
    cmd.extend(str(x) for x in (extra or []))
    _run(cmd)
    return out


def load_mag_taxonomy(path: Path) -> Dict[str, str]:
    out: Dict[str, str] = {}
    if not path.is_file():
        return out
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for rec in reader:
            mag = rec.get("mag_id") or rec.get("user_genome") or rec.get("Name") or ""
            mag = str(mag).strip()
            taxid = rec.get("taxid") or rec.get("taxID") or "0"
            if mag:
                tid = str(taxid).split(".")[0] or "0"
                prev = out.get(mag)
                if prev and prev != "0" and tid == "0":
                    continue
                out[mag] = tid
    return out


def rewrite_mag_taxonomy_ncbi(path: Any) -> Path:
    """Fill NCBI taxids from GTDB-Tk lineage strings already on disk."""
    src = as_path(path)
    if not src.is_file():
        raise FileNotFoundError(src)
    rows: List[Dict[str, str]] = []
    with src.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = list(reader.fieldnames or ["mag_id", "taxid", "lineage"])
        if "taxid" not in fieldnames:
            fieldnames.insert(1, "taxid")
        for rec in reader:
            mag = rec.get("mag_id") or rec.get("user_genome") or rec.get("Name") or ""
            mag = str(mag).strip()
            if not mag:
                continue
            lineage = rec.get("lineage") or rec.get("classification") or ""
            rec["mag_id"] = mag
            rec["lineage"] = str(lineage)
            rec["taxid"] = _gtdb_taxid_from_lineage(str(lineage))
            if rec["taxid"] == "0" and not str(lineage).strip():
                continue
            rows.append(rec)
    with src.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    return src


def run_taxon_quantifier(mag_abundance: str, taxonomy: str, dest: str) -> Path:
    tax = load_mag_taxonomy(as_path(taxonomy))
    counts: Dict[str, int] = defaultdict(int)
    abund = as_path(mag_abundance)
    if abund.is_file():
        with abund.open(newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            field = "N"
            if reader.fieldnames:
                for name in reader.fieldnames:
                    if name == "N" or str(name).startswith("N_"):
                        field = name
                        break
            for rec in reader:
                mag = rec.get("mag_id") or rec.get("Genome") or rec.get("genome") or ""
                n = rec.get(field) or "0"
                try:
                    val = int(float(n))
                except ValueError:
                    val = 0
                counts[tax.get(str(mag), "0")] += val
    out = as_path(dest)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = ["taxid\tN_1"]
    if not counts:
        rows.append("0\t0")
    else:
        for taxid, n in sorted(counts.items()):
            rows.append(f"{taxid}\t{n}")
    out.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return out


def _alignment_needs_samtools(path: Path) -> bool:
    """True for BAM/CRAM. BGZF BAM starts with gzip magic, not ``BAM\\x01``.

    A ``.bam`` suffix is not enough: tests and some pipelines write SAM text
    under that name, and sending it to ``samtools view`` fails.
    """
    if not path.is_file():
        return False
    head = path.read_bytes()[:4]
    if head == b"BAM\x01" or head[:4] == b"CRAM" or head[:2] == b"\x1f\x8b":
        return True
    return path.suffix.lower() == ".cram"


def _iter_sam_alignments(bam: Path) -> Iterable[Tuple[str, str]]:
    if not bam.is_file():
        return
    if _alignment_needs_samtools(bam):
        view = subprocess.check_output(
            [which_tool("samtools"), "view", "-F", "2308", str(bam)],
            text=True,
            errors="replace",
        )
        lines = view.splitlines()
    else:
        lines = bam.read_text(encoding="utf-8", errors="replace").splitlines()
    for line in lines:
        if not line or line.startswith("@"):
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        yield parts[0], parts[2]


def run_read_assigner(bam: str, taxonomy: str, dest: str, mag_dir: str = "") -> Path:
    tax = load_mag_taxonomy(as_path(taxonomy))
    mapping = contig_to_mag_map(as_path(mag_dir)) if mag_dir else {}
    dest_path = as_path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for seq, rname in _iter_sam_alignments(as_path(bam)):
        if rname in {"*", ""}:
            continue
        mag = mapping.get(rname, rname.split("~", 1)[0])
        taxid = tax.get(mag, tax.get(rname, "0"))
        rows.append(f"{seq}\t{taxid}")
    dest_path.write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")
    return dest_path


def run_assembly_annotator(
    r1: str,
    r2: str,
    dest: str,
    db: str = "",
    threads: int = 1,
    *,
    assembler: str = "megahit",
    gene_caller: str = "prodigal",
    binners: Optional[Sequence[str]] = None,
    binner_qc: str = "checkm2",
    binner_combine: str = "dastool",
    mag_taxonomy: str = "gtdbtk",
    aligner: str = "minimap2",
    mag_quantifier: str = "coverm",
    work_dir: str = "",
) -> Path:
    """Full nested pipeline; writes annotator TSV at dest."""
    out = as_path(dest)
    work = as_path(work_dir) if work_dir else out.with_suffix("")
    if str(work) == str(out):
        work = out.parent / (out.name + ".assembly")
    work.mkdir(parents=True, exist_ok=True)
    slots = parse_slot_extra(os.environ.get("SAMOVAR_ASSEMBLY_EXTRA", ""))
    if slots.get("assembler"):
        assembler = str(slots["assembler"])
    if slots.get("gene_caller"):
        gene_caller = str(slots["gene_caller"])
    names = list(binners or ["anvio"])
    if slots.get("binners"):
        names = list(slots["binners"])
    if slots.get("binner_qc"):
        binner_qc = str(slots["binner_qc"])
    if slots.get("binner_combine"):
        binner_combine = str(slots["binner_combine"])
    if slots.get("mag_taxonomy"):
        mag_taxonomy = str(slots["mag_taxonomy"])
    if slots.get("aligner"):
        aligner = str(slots["aligner"])
    if slots.get("mag_quantifier"):
        mag_quantifier = str(slots["mag_quantifier"])
    contigs = work / "contigs.fa"
    run_assembler(r1, r2, str(contigs), threads=threads, name=assembler)
    genes = work / "genes"
    try:
        run_gene_caller(str(contigs), str(genes), threads=threads, name=gene_caller)
    except (FileNotFoundError, subprocess.CalledProcessError, OSError):
        genes.mkdir(parents=True, exist_ok=True)
    bin_dirs = []
    for name in names:
        folder = work / "bins" / name
        try:
            run_binner(str(contigs), str(folder), r1=r1, r2=r2, threads=threads, name=name)
        except (FileNotFoundError, subprocess.CalledProcessError, OSError):
            _write_identity_mag(contigs, folder)
        bin_dirs.append(str(folder))
        qc_path = work / "bins_qc" / f"{name}.tsv"
        try:
            run_binner_qc(str(folder), str(qc_path), db="", threads=threads, name=binner_qc)
        except (FileNotFoundError, subprocess.CalledProcessError):
            dummy_qc = as_path(qc_path)
            dummy_qc.parent.mkdir(parents=True, exist_ok=True)
            dummy_qc.write_text("mag_id\tcompleteness\tcontamination\tscore\n", encoding="utf-8")
    selected = work / "bins_selected"
    try:
        run_binner_combine(
            bin_dirs,
            str(selected),
            contigs=str(contigs),
            threads=threads,
            name=binner_combine,
        )
    except (FileNotFoundError, subprocess.CalledProcessError):
        run_binner_combine(bin_dirs, str(selected), name="identity")
    mag_genes = work / "mag_genes"
    try:
        run_gene_caller(str(selected), str(mag_genes), threads=threads, name=gene_caller)
    except (FileNotFoundError, subprocess.CalledProcessError, OSError):
        mag_genes.mkdir(parents=True, exist_ok=True)
    tax_table = work / "mag_taxonomy.tsv"
    try:
        run_mag_taxonomy(str(selected), str(tax_table), db=db, threads=threads, name=mag_taxonomy)
    except (FileNotFoundError, subprocess.CalledProcessError):
        rows = ["mag_id\ttaxid\tlineage"]
        for fasta in _iter_mag_fastas(selected):
            rows.append(f"{fasta.stem}\t0\t")
        tax_table.write_text("\n".join(rows) + "\n", encoding="utf-8")
    dest_path = as_path(out)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    if not any(_fasta_has_records(p) for p in _iter_mag_fastas(selected)):
        dest_path.write_text("", encoding="utf-8")
        return dest_path
    bam = work / "reads_to_mag.bam"
    run_aligner(r1, r2, str(selected), str(bam), threads=threads, name=aligner)
    mag_abund = work / "mag_abundance.tsv"
    run_mag_quantifier(str(bam), str(selected), str(mag_abund), name=mag_quantifier)
    taxon_abund = work / "taxon_abundance.tsv"
    run_taxon_quantifier(str(mag_abund), str(tax_table), str(taxon_abund))
    return run_read_assigner(str(bam), str(tax_table), str(out), mag_dir=str(selected))


def _stage_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="python -m samovar.assembly_profiling")
    parser.add_argument("stage", nargs="?", default="annotator", choices=STAGE_NAMES)
    parser.add_argument("-i", dest="r1", default="")
    parser.add_argument("-I", dest="r2", default="")
    parser.add_argument("-c", dest="contigs", default="")
    parser.add_argument("-r", dest="reference", default="")
    parser.add_argument("-b", dest="bam", default="")
    parser.add_argument("-a", dest="abundance", default="")
    parser.add_argument("-x", dest="taxonomy", default="")
    parser.add_argument("-q", dest="qc", default="")
    parser.add_argument("-d", dest="db", default="")
    parser.add_argument("-o", dest="output", default="")
    parser.add_argument("-t", dest="threads", default="1")
    parser.add_argument("--assembler", default="megahit")
    parser.add_argument("--gene-caller", dest="gene_caller", default="prodigal")
    parser.add_argument("--binner", action="append", dest="binners")
    parser.add_argument("--binner-qc", dest="binner_qc", default="checkm2")
    parser.add_argument("--binner-combine", dest="binner_combine", default="dastool")
    parser.add_argument("--mag-taxonomy", dest="mag_taxonomy", default="gtdbtk")
    parser.add_argument("--aligner", default="minimap2")
    parser.add_argument("--mag-quantifier", dest="mag_quantifier", default="coverm")
    parser.add_argument("--work-dir", dest="work_dir", default="")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _stage_parser().parse_args(list(argv) if argv is not None else None)
    threads = int(args.threads or 1)
    stage = args.stage
    if stage == "assembler":
        run_assembler(args.r1, args.r2, args.output, threads=threads, name=args.assembler)
    elif stage == "gene_caller":
        run_gene_caller(
            args.contigs or args.r1,
            args.output,
            threads=threads,
            name=args.gene_caller,
        )
    elif stage == "binner":
        names = args.binners or ["anvio"]
        run_binner(args.contigs, args.output, r1=args.r1, r2=args.r2, threads=threads, name=names[0])
    elif stage == "binner_qc":
        run_binner_qc(args.contigs, args.output, db=args.db, threads=threads, name=args.binner_qc)
    elif stage == "binner_combine":
        dirs = [p for p in str(args.contigs).split(",") if p]
        run_binner_combine(dirs, args.output, qc_table=args.qc, threads=threads, name=args.binner_combine)
    elif stage == "mag_taxonomy":
        run_mag_taxonomy(args.contigs, args.output, db=args.db, threads=threads, name=args.mag_taxonomy)
    elif stage == "aligner":
        run_aligner(args.r1, args.r2, args.reference, args.output, threads=threads, name=args.aligner)
    elif stage == "mag_quantifier":
        run_mag_quantifier(args.bam, args.reference, args.output, name=args.mag_quantifier)
    elif stage == "taxon_quantifier":
        run_taxon_quantifier(args.abundance, args.taxonomy, args.output)
    elif stage == "read_assigner":
        run_read_assigner(args.bam, args.taxonomy, args.output, mag_dir=args.reference)
    else:
        run_assembly_annotator(
            args.r1,
            args.r2,
            args.output,
            db=args.db,
            threads=threads,
            assembler=args.assembler,
            gene_caller=args.gene_caller,
            binners=args.binners,
            binner_qc=args.binner_qc,
            binner_combine=args.binner_combine,
            mag_taxonomy=args.mag_taxonomy,
            aligner=args.aligner,
            mag_quantifier=args.mag_quantifier,
            work_dir=args.work_dir,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
