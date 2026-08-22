"""
Module for converting tables to simulated reads using ISS-like functionality.
"""

import os
import re
import glob
import gzip
import random
import pandas as pd
import subprocess
from .genome_fetcher import fetch_genome, default_entrez_email
from typing import Dict, Iterable, List, Optional, Sequence, Tuple, Union
import yaml
import json
import tempfile
import warnings
import shutil
from samovar.paths import annotation_regenerate_r, iss_executable, load_config, resolve_executable


SKIP_TAXIDS = {"0", "nan", "None", ""}
CONTIG_SPACER = "N" * 500


def _iss_cmd(*args: str) -> List[str]:
    return [iss_executable(), *args]


def _run_iss(cmd: Sequence[str]) -> None:
    result = subprocess.run(list(cmd), check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"ISS command failed (exit {result.returncode}): {' '.join(cmd)}"
        )


def parse_annotation_table(table_path: str) -> pd.DataFrame:
    """
    Parse annotation table and return DataFrame with taxid and counts for each annotation method.
    
    Args:
        table_path: Path to taxonomy table
        
    Returns:
        DataFrame with columns: taxid, N_k1, N_k2, ... where each N_k column contains
        the count of occurrences for that annotation method
    """
    # Read the input table
    df = pd.read_csv(table_path, sep=",")
    
    # Get all taxid columns (they contain 'taxid' in their name)
    taxid_cols = [col for col in df.columns if 'taxid' in col.lower()]
    
    # Create a new DataFrame to store results
    result = pd.DataFrame()
    
    # For each taxid column, calculate counts
    for col in taxid_cols:
        # Get the annotation name from the column name (e.g., 'k1' from 'taxID_k1')
        ann_name = col.split('_')[-2]
        
        # Count occurrences of each taxid
        taxid_counts = df.groupby(col).size().reset_index()
        taxid_counts.columns = ["taxid", f"N_{ann_name}"]
        
        # Merge with result DataFrame
        if result.empty:
            result = taxid_counts
        else:
            result = result.merge(taxid_counts, on='taxid', how='outer')
    
    # Fill NaN values with 0
    result = result.fillna(0)

    if result.empty or "taxid" not in result.columns:
        return pd.DataFrame(columns=["taxid"])

    # Convert taxid to string
    result['taxid'] = result['taxid'].astype(str)
    
    return result

def get_genome_file(genome_dir: str, taxid: str) -> str:
    """
    Get the path to a genome file, checking multiple possible extensions.
    
    Args:
        genome_dir: Directory containing genome files
        taxid: Taxonomy ID of the genome
        
    Returns:
        Path to the genome file if found, None otherwise
    """
    extensions = ['-processed.fa', '-processed.fna', '-processed.fasta',
                  '.fa', '.fna', '.fasta', '.fa.gz', '.fna.gz', '.fasta.gz']
    
    for ext in extensions:
        genome_file = os.path.join(genome_dir, f"{taxid}{ext}")
        if os.path.exists(genome_file):
            return genome_file
    return None


def _n_columns(df: pd.DataFrame) -> List[str]:
    return [col for col in df.columns if col != "taxid" and "n" in col.lower()]


def _annotator_from_n_col(column: str) -> str:
    match = re.search(r"N_(.*?)(?:_[0-9]*)?$", column)
    return match.group(1) if match else "any"


def _safe_fasta_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9.-]", "_", str(value))
    return cleaned or "genome"


def _genome_id_from_path(path: str) -> str:
    name = os.path.basename(path)
    for ext in (
        ".fasta.gz", ".fna.gz", ".fa.gz",
        "-processed.fasta", "-processed.fna", "-processed.fa",
        ".fasta", ".fna", ".fa",
    ):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return _safe_fasta_id(name)


def _open_text(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def _read_fasta_concat(path: str) -> str:
    seqs = []
    current = []
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if current:
                    seqs.append("".join(current))
                    current = []
            else:
                current.append(line.strip())
        if current:
            seqs.append("".join(current))
    return CONTIG_SPACER.join(seqs)


def _scale_amounts(amount: Sequence[int], total_amount: Optional[int]) -> List[int]:
    scaled = [max(0, int(n)) for n in amount]
    total = sum(scaled)
    if total_amount is not None and total > 0:
        scaled = [int(n / total * total_amount) for n in scaled]
    return scaled


def _write_empty_fastq_pair(r1: str, r2: str) -> None:
    for path in (r1, r2):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open(path, "w") as handle:
            handle.write("\n")


def _cleanup_iss_tmp(prefix: str) -> None:
    directory = os.path.dirname(prefix) or "."
    base = os.path.basename(prefix)
    for path in glob.glob(os.path.join(directory, f"{base}*")):
        name = os.path.basename(path)
        if name.endswith("_R1.fastq") or name.endswith("_R2.fastq"):
            continue
        if os.path.isfile(path):
            try:
                os.remove(path)
            except OSError:
                pass


def _iss_record_id(genome_id: str) -> str:
    """FASTA/read ID that ISS can count and that true-label parsing can recover."""
    gid = str(genome_id)
    if re.fullmatch(r"[0-9]+", gid):
        return f"taxid:{gid}"
    return _safe_fasta_id(gid)


def _source_id_from_header(header: str) -> str:
    rid = header.lstrip("@").split()[0].split("/")[0]
    match = re.search(r"taxid:([0-9]+)", rid)
    if match:
        return match.group(1)
    return rid.split("_")[0]


def _iter_fastq_pairs(r1_path: str, r2_path: str):
    with open(r1_path) as handle_r1, open(r2_path) as handle_r2:
        while True:
            rec1 = [handle_r1.readline() for _ in range(4)]
            rec2 = [handle_r2.readline() for _ in range(4)]
            if not rec1[0] or not rec2[0]:
                break
            if not rec1[0].strip():
                continue
            yield rec1, rec2


def _write_fastq_records(path: str, records: Iterable[List[str]]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    wrote = False
    with open(path, "w") as handle:
        for record in records:
            handle.write("".join(record))
            wrote = True
        if not wrote:
            handle.write("\n")


def _allocate_counts(available: int, requested: Dict[str, int]) -> Dict[str, int]:
    samples = list(requested.keys())
    total_requested = sum(max(0, int(requested[s])) for s in samples)
    allocated = {s: 0 for s in samples}
    if available <= 0 or not samples:
        return allocated
    if total_requested <= 0:
        share = available // len(samples)
        remainder = available - share * len(samples)
        for i, sample in enumerate(samples):
            allocated[sample] = share + (1 if i < remainder else 0)
        return allocated

    remaining = available
    for i, sample in enumerate(samples):
        if i == len(samples) - 1:
            allocated[sample] = remaining
        else:
            n = int(round(max(0, int(requested[sample])) / total_requested * available))
            n = min(max(0, n), remaining)
            allocated[sample] = n
            remaining -= n
    return allocated


def write_combined_genomes(
    genome_files: Sequence[str],
    dest_fasta: str,
    genome_ids: Optional[Sequence[str]] = None,
) -> List[str]:
    """Write one FASTA record per genome with a stable ID for ISS + later splitting."""
    os.makedirs(os.path.dirname(dest_fasta) or ".", exist_ok=True)
    if genome_ids is None:
        genome_ids = [_genome_id_from_path(path) for path in genome_files]

    written_ids = []
    with open(dest_fasta, "w") as out:
        for genome_file, genome_id in zip(genome_files, genome_ids):
            if genome_file is None or not os.path.exists(genome_file):
                continue
            seq = _read_fasta_concat(genome_file)
            if not seq:
                continue
            fasta_id = _iss_record_id(genome_id)
            out.write(f">{fasta_id}\n{seq}\n")
            written_ids.append(fasta_id)
    return written_ids


def write_readcount_file(
    dest_path: str,
    genome_ids: Sequence[str],
    amount: Sequence[int],
) -> int:
    """Write an ISS --readcount_file (record_id<TAB>n_reads). Returns total reads."""
    os.makedirs(os.path.dirname(dest_path) or ".", exist_ok=True)
    total = 0
    with open(dest_path, "w") as handle:
        for genome_id, n_reads in zip(genome_ids, amount):
            n_reads = max(0, int(n_reads))
            if n_reads <= 0:
                continue
            handle.write(f"{genome_id}\t{n_reads}\n")
            total += n_reads
    return total


def split_paired_fastq_by_counts(
    r1_path: str,
    r2_path: str,
    sample_counts: Dict[str, int],
    output_paths: Dict[str, Tuple[str, str]],
) -> Dict[str, int]:
    """
    Split a paired FASTQ into samples using sequential chunks proportional to counts.
    """
    if not os.path.exists(r1_path) or not os.path.exists(r2_path):
        for sample, (r1_out, r2_out) in output_paths.items():
            _write_empty_fastq_pair(r1_out, r2_out)
        return {sample: 0 for sample in output_paths}

    pairs = list(_iter_fastq_pairs(r1_path, r2_path))
    allocated = _allocate_counts(len(pairs), sample_counts)
    offset = 0
    written = {}
    for sample in sample_counts:
        n = allocated.get(sample, 0)
        chunk = pairs[offset: offset + n]
        offset += n
        r1_out, r2_out = output_paths[sample]
        _write_fastq_records(r1_out, (rec[0] for rec in chunk))
        _write_fastq_records(r2_out, (rec[1] for rec in chunk))
        written[sample] = len(chunk)
    return written


def split_metagenome_to_samples(
    r1_path: str,
    r2_path: str,
    sample_source_counts: Dict[str, Dict[str, int]],
    output_dir: str,
    annotator_name: str,
) -> Dict[str, int]:
    """
    Split a mixed metagenome FASTQ into per-sample files, preserving source composition.

    Read headers are expected to start with the genome/taxid ID used during generation
    (`{source_id}_{n}_{cpu}/1`).
    """
    samples = list(sample_source_counts.keys())
    output_handles = {}
    written = {sample: 0 for sample in samples}

    for sample in samples:
        r1_out = os.path.join(output_dir, f"{sample}_{annotator_name}_R1.fastq")
        r2_out = os.path.join(output_dir, f"{sample}_{annotator_name}_R2.fastq")
        os.makedirs(output_dir, exist_ok=True)
        output_handles[sample] = (open(r1_out, "w"), open(r2_out, "w"))

    try:
        if not os.path.exists(r1_path) or not os.path.exists(r2_path):
            return written

        buckets: Dict[str, List[Tuple[List[str], List[str]]]] = {}
        for rec1, rec2 in _iter_fastq_pairs(r1_path, r2_path):
            source_id = _source_id_from_header(rec1[0])
            buckets.setdefault(source_id, []).append((rec1, rec2))

        all_sources = set()
        for counts in sample_source_counts.values():
            all_sources.update(counts.keys())
        all_sources.update(buckets.keys())

        for source_id in all_sources:
            available = buckets.get(source_id, [])
            requested = {
                sample: int(sample_source_counts.get(sample, {}).get(source_id, 0))
                for sample in samples
            }
            allocated = _allocate_counts(len(available), requested)
            offset = 0
            for sample in samples:
                n = allocated.get(sample, 0)
                chunk = available[offset: offset + n]
                offset += n
                r1_handle, r2_handle = output_handles[sample]
                for rec1, rec2 in chunk:
                    r1_handle.write("".join(rec1))
                    r2_handle.write("".join(rec2))
                    written[sample] += 1
    finally:
        for r1_handle, r2_handle in output_handles.values():
            for handle in (r1_handle, r2_handle):
                if handle.tell() == 0:
                    handle.write("\n")
                handle.close()

    return written


def generate_reads_genome(
    genome_file: str,
    output_file: str,
    amount: int,
    read_length: int = 150,
    model: str = "hiseq",
    cpus: int = 2,
    seed: Optional[int] = None,
) -> None:
    """
    Generate simulated reads from a genome file.
    
    Args:
        genome_file: Path to genome FASTA file
        output_file: Path to output FASTQ file
        amount: Number of reads to generate
        read_length: Length of reads to generate
        model: Model to use for simulation
        cpus: ISS worker count
        seed: Optional ISS seed
    """
    
    if genome_file is not None and os.path.exists(genome_file):
        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        cmd = _iss_cmd(
            "generate",
            "--genomes", genome_file,
            "--model", model,
            "--output", output_file,
            "--n_reads", str(int(amount)),
            "--cpus", str(cpus),
        )
        if seed is not None:
            cmd.extend(["--seed", str(seed)])
        _run_iss(cmd)
        _cleanup_iss_tmp(output_file)


def generate_reads_metagenome(
    genome_files: List[str],
    output_dir: str,
    amount: List[int],
    read_length: int = 150,
    total_amount: int = None,
    sample_name: str = "merged",
    annotator_name: str = None,
    model: str = "hiseq",
    genome_ids: Optional[Sequence[str]] = None,
    cpus: int = 2,
    seed: Optional[int] = None,
) -> Tuple[str, str]:
    """
    Generate simulated reads from a metagenome in a single ISS invocation.

    All genomes are combined into one community (readcount-weighted), then ISS
    writes the mixed R1/R2 FASTQ directly — no per-genome ISS calls.
    """
    amount = _scale_amounts(amount, total_amount)
    if annotator_name is None:
        annotator_name = "any"
    if genome_ids is None:
        genome_ids = [_genome_id_from_path(path) for path in genome_files]

    os.makedirs(output_dir, exist_ok=True)
    output_prefix = os.path.join(output_dir, f"{sample_name}_{annotator_name}")
    output_file_R1 = f"{output_prefix}_R1.fastq"
    output_file_R2 = f"{output_prefix}_R2.fastq"

    pairs = [
        (path, gid, n)
        for path, gid, n in zip(genome_files, genome_ids, amount)
        if path is not None and os.path.exists(path) and int(n) > 0
    ]
    if not pairs:
        _write_empty_fastq_pair(output_file_R1, output_file_R2)
        return output_file_R1, output_file_R2

    work_dir = os.path.join(output_dir, ".iss_full")
    os.makedirs(work_dir, exist_ok=True)
    combined_fasta = os.path.join(work_dir, f"{sample_name}_{annotator_name}.genomes.fasta")
    readcount_path = os.path.join(work_dir, f"{sample_name}_{annotator_name}.readcount.txt")

    written_ids = write_combined_genomes(
        [p[0] for p in pairs],
        combined_fasta,
        [p[1] for p in pairs],
    )
    amount_by_id = {_iss_record_id(gid): n_reads for _, gid, n_reads in pairs}
    total_reads = write_readcount_file(
        readcount_path,
        written_ids,
        [amount_by_id[gid] for gid in written_ids],
    )
    if total_reads <= 0 or os.path.getsize(combined_fasta) == 0:
        _write_empty_fastq_pair(output_file_R1, output_file_R2)
        return output_file_R1, output_file_R2

    cmd = _iss_cmd(
        "generate",
        "--genomes", combined_fasta,
        "--readcount_file", readcount_path,
        "--model", model,
        "--output", output_prefix,
        "--cpus", str(cpus),
    )
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    _run_iss(cmd)
    _cleanup_iss_tmp(output_prefix)

    if not os.path.exists(output_file_R1) or not os.path.exists(output_file_R2):
        raise RuntimeError(
            f"ISS did not produce reads for {output_prefix}"
        )

    return output_file_R1, output_file_R2

def regenerate_metagenome(
    genome_files: List[str],
    output_dir: str,
    amount: List[int],
    read_length: int = 150,
    total_amount: int = None,
    sample_name: str = None,
    annotator_name: str = None,
    model: str = "hiseq",
    genome_ids: Optional[Sequence[str]] = None,
    cpus: int = 2,
    seed: Optional[int] = None,
) -> Tuple[str, str]:
    """
    Regenerate metagenome reads from a list of genome files.
    """
    return generate_reads_metagenome(
        genome_files=genome_files,
        output_dir=output_dir,
        amount=amount,
        read_length=read_length,
        total_amount=total_amount,
        model=model,
        sample_name=sample_name,
        annotator_name=annotator_name,
        genome_ids=genome_ids,
        cpus=cpus,
        seed=seed,
    )

def process_annotation_table(
    table_path: str,
    genome_dir: str,
    output_dir: str,
    total_amount: int = None,
    email: Optional[str] = None,
    reference_only: bool = True,
    model: str = "hiseq",
    read_length: int = 150,
    sample_name: str = None
) -> None:
    """
    Process taxonomy table and generate simulated reads for each taxid. Inherits from process_abundance_table.
    
    Args:
        table_path: Path to taxonomy table file
        genome_dir: Directory containing genome files
        output_dir: Directory to write output files
        read_length: Length of reads to generate
        email: Email for NCBI Entrez
        reference_only: If True, only fetch reference genome
        model: Model to use for simulation with ISS
        sample_name: Name of the sample
        total_amount: Total number of reads to generate
    """
    email = email or default_entrez_email()
    abundance_table = parse_annotation_table(table_path)
    process_abundance_table(
        table=abundance_table,
        genome_dir=genome_dir,
        output_dir=output_dir,
        total_amount=total_amount,
        email=email,
        reference_only=reference_only,
        model=model,
        read_length=read_length,
        sample_name=sample_name
    )
        
def _resolve_genomes_for_taxids(
    taxids: Iterable[str],
    genome_dir: str,
    email: str,
    reference_only: bool,
    max_genomes: Optional[int] = None,
) -> Dict[str, str]:
    available = {}
    for taxid in taxids:
        if max_genomes is not None and len(available) >= max_genomes:
            break
        taxid = str(taxid).split(".")[0]
        if taxid in SKIP_TAXIDS or taxid in available:
            continue
        genome_file = get_genome_file(genome_dir, taxid)
        if genome_file is None:
            try:
                genome_file = fetch_genome(
                    taxid, genome_dir, email, reference_only=reference_only
                )
            except Exception as exc:
                warnings.warn(f"Failed to fetch genome for taxid {taxid}: {exc}")
                genome_file = None
        if genome_file is not None:
            available[taxid] = genome_file
    return available


def _emit_empty_for_annotators(
    output_dir: str,
    sample_names: Sequence[str],
    annotator_names: Sequence[str],
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    for sample_name in sample_names:
        for annotator_name in annotator_names:
            _write_empty_fastq_pair(
                os.path.join(output_dir, f"{sample_name}_{annotator_name}_R1.fastq"),
                os.path.join(output_dir, f"{sample_name}_{annotator_name}_R2.fastq"),
            )


def process_abundance_table(
    table: Union[str, pd.DataFrame],
    genome_dir: str,
    output_dir: str,
    total_amount: int = None,
    email: Optional[str] = None,
    reference_only: bool = True,
    model: str = "hiseq",
    read_length: int = 150,
    sample_name: str = None
) -> pd.DataFrame:
    """
    Process abundance table and generate simulated reads for each taxid.
    
    Args:
        table: Path to taxonomy table file or pandas DataFrame
        genome_dir: Directory containing genome files
        output_dir: Directory to write output files
        read_length: Length of reads to generate
        email: Email for NCBI Entrez
        reference_only: If True, only fetch reference genome
        model: Model to use for simulation with ISS
        sample_name: Name of the sample
        total_amount: Total number of reads to generate
        
    Returns:
        DataFrame containing the filtered abundance table with only available genomes
    """
    if isinstance(table, str):
        abundance_table = pd.read_csv(table, sep=",")
    else:
        abundance_table = table
    
    if sample_name is None:
        if isinstance(table, str):
            sample_name = os.path.basename(table).split(".")[0]
        else:
            sample_name = "merged"

    email = email or default_entrez_email()

    if "taxid" not in abundance_table.columns:
        _emit_empty_for_annotators(output_dir, [sample_name], ["any"])
        warnings.warn("Abundance table has no taxid column; emitted empty FASTQ files")
        return None

    abundance_table = abundance_table.copy()
    abundance_table["taxid"] = abundance_table["taxid"].astype(str).str.split(".").str[0]
    available_genomes = _resolve_genomes_for_taxids(
        abundance_table["taxid"], genome_dir, email, reference_only
    )

    N_cols = _n_columns(abundance_table)
    annotators = [_annotator_from_n_col(col) for col in N_cols] or ["any"]

    # Gracefully handle all-unclassified/unknown samples.
    if not available_genomes:
        _emit_empty_for_annotators(output_dir, [sample_name], annotators)
        warnings.warn("No genome files available for any taxid; emitted empty FASTQ files")
        return None

    filtered_table = abundance_table[abundance_table["taxid"].isin(available_genomes)]
    N_cols = _n_columns(filtered_table)
    if N_cols:
        filtered_table = filtered_table.groupby("taxid", as_index=False)[N_cols].sum()

    os.makedirs(output_dir, exist_ok=True)
    for N_annotator in N_cols:
        annotator_name = _annotator_from_n_col(N_annotator)
        genome_ids = filtered_table["taxid"].tolist()
        genome_files = [available_genomes[taxid] for taxid in genome_ids]
        amount = filtered_table[N_annotator].tolist()
        regenerate_metagenome(
            genome_files=genome_files,
            output_dir=output_dir,
            amount=amount,
            read_length=read_length,
            total_amount=total_amount,
            sample_name=sample_name,
            model=model,
            annotator_name=annotator_name,
            genome_ids=genome_ids,
        )

    return filtered_table


def _sample_tables_from_abundance_dir(
    abundance_dir: str,
    sample_names_hint: Optional[Sequence[str]] = None,
) -> Dict[str, pd.DataFrame]:
    """Build per-sample abundance tables from regenerated annotator CSVs."""
    from pathlib import Path

    from samovar.regenerate import sample_names_from_abundance_columns

    by_annotator: Dict[str, pd.DataFrame] = {}
    for path in sorted(Path(abundance_dir).glob("*.csv")):
        by_annotator[path.stem] = pd.read_csv(path)
    if not by_annotator:
        return {}

    ref = next(iter(by_annotator.values()))
    n_cols = _n_columns(ref)
    if not n_cols:
        return {}

    sample_names = sample_names_from_abundance_columns(
        n_cols,
        list(sample_names_hint) if sample_names_hint else None,
    )
    sample_tables: Dict[str, pd.DataFrame] = {}
    for idx, n_col in enumerate(n_cols):
        sample_name = sample_names[idx]
        parts: List[pd.DataFrame] = []
        for ann_name, df in by_annotator.items():
            if n_col not in df.columns:
                continue
            sub = df[["taxid", n_col]].copy()
            sub.columns = ["taxid", f"N_{ann_name}"]
            parts.append(sub)
        if not parts:
            continue
        merged = parts[0]
        for part in parts[1:]:
            merged = merged.merge(part, on="taxid", how="outer")
        merged = merged.fillna(0)
        merged["taxid"] = merged["taxid"].astype(str)
        sample_tables[sample_name] = merged
    return sample_tables


def process_annotation_tables(
    table_paths: Sequence[str],
    genome_dir: str,
    output_dir: str,
    total_amount: int = None,
    email: Optional[str] = None,
    reference_only: bool = True,
    model: str = "hiseq",
    read_length: int = 150,
    sample_names: Optional[Sequence[str]] = None,
    cpus: int = 2,
    seed: Optional[int] = None,
    max_genomes: Optional[int] = None,
    annotation_dir: Optional[str] = None,
    regeneration_config: Optional[dict] = None,
) -> None:
    """
    Generate one full metagenome per annotator, then split reads into samples.

    This is much faster than calling ISS once per sample (and per genome).

    Args:
        max_genomes: If set, only resolve/fetch the top-N taxids by total
            abundance across samples (critical for large public DBs).
        annotation_dir: Directory of per-sample ``*.annotation.csv`` files.
            Required when ``regeneration_config`` selects a generative mode.
        regeneration_config: SamovaR-style regeneration settings (``regeneration_mode``,
            ``N``, ``N_reads``, ``seed``, etc.).
    """
    from samovar.regenerate import normalize_regeneration_mode, write_samovar_config_defaults

    regen_cfg = write_samovar_config_defaults(dict(regeneration_config or {}))
    mode = normalize_regeneration_mode(regen_cfg.get("regeneration_mode"))

    table_paths = list(table_paths)
    email = email or default_entrez_email()
    if sample_names is None:
        sample_names = [os.path.basename(path).split(".")[0] for path in table_paths]
    else:
        sample_names = [
            name or os.path.basename(path).split(".")[0]
            for name, path in zip(sample_names, table_paths)
        ]

    sample_tables: Dict[str, pd.DataFrame] = {}
    if mode != "preserve" and annotation_dir:
        regen_dir = os.path.join(output_dir, ".regenerated_abundance")
        cfg_out = dict(regen_cfg)
        cfg_out["output_dir"] = regen_dir
        cfg_out["regeneration_mode"] = mode
        with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as tmp:
            yaml.dump(cfg_out, tmp)
            tmp_config = tmp.name
        try:
            samovar_annotation_regenerate(
                annotation_dir=str(annotation_dir),
                config_samovar=tmp_config,
                output_dir=regen_dir,
            )
        finally:
            os.unlink(tmp_config)
        sample_tables = _sample_tables_from_abundance_dir(regen_dir, sample_names)
        sample_names = list(sample_tables.keys())
    else:
        for path, sample_name in zip(table_paths, sample_names):
            table = parse_annotation_table(path)
            if not table.empty and "taxid" in table.columns:
                table = table.copy()
                table["taxid"] = table["taxid"].astype(str).str.split(".").str[0]
            sample_tables[sample_name] = table

    if not sample_tables:
        _emit_empty_for_annotators(output_dir, sample_names or ["1"], ["any"])
        warnings.warn("No annotation or abundance tables to process; emitted empty FASTQ files")
        return

    annotator_cols: Dict[str, str] = {}
    taxid_totals: Dict[str, int] = {}
    for table in sample_tables.values():
        for col in _n_columns(table):
            annotator_cols[_annotator_from_n_col(col)] = col
            if not table.empty and "taxid" in table.columns:
                for taxid, n_reads in zip(table["taxid"].astype(str), table[col]):
                    taxid_totals[taxid] = taxid_totals.get(taxid, 0) + int(n_reads)

    # Prefer abundant taxa when capping genome downloads
    ranked_taxids = sorted(taxid_totals, key=taxid_totals.get, reverse=True)

    annotators = list(annotator_cols.keys()) or ["any"]
    available_genomes = _resolve_genomes_for_taxids(
        ranked_taxids,
        genome_dir,
        email,
        reference_only,
        max_genomes=max_genomes,
    )
    if not available_genomes:
        _emit_empty_for_annotators(output_dir, sample_names, annotators)
        warnings.warn("No genome files available for any taxid; emitted empty FASTQ files")
        return

    os.makedirs(output_dir, exist_ok=True)
    pool_dir = os.path.join(output_dir, ".iss_full")
    os.makedirs(pool_dir, exist_ok=True)

    for annotator_name in annotators:
        sample_source_counts: Dict[str, Dict[str, int]] = {}
        total_by_taxid: Dict[str, int] = {}
        for sample_name, table in sample_tables.items():
            counts: Dict[str, int] = {}
            n_col = next(
                (col for col in _n_columns(table) if _annotator_from_n_col(col) == annotator_name),
                None,
            )
            if n_col is not None and not table.empty:
                amounts = table[n_col].tolist()
                taxids = table["taxid"].tolist()
                amounts = _scale_amounts(amounts, total_amount)
                for taxid, n_reads in zip(taxids, amounts):
                    taxid = str(taxid)
                    if taxid not in available_genomes or n_reads <= 0:
                        continue
                    counts[taxid] = counts.get(taxid, 0) + int(n_reads)
                    total_by_taxid[taxid] = total_by_taxid.get(taxid, 0) + int(n_reads)
            sample_source_counts[sample_name] = counts

        genome_ids = [taxid for taxid, n_reads in total_by_taxid.items() if n_reads > 0]
        if not genome_ids:
            _emit_empty_for_annotators(output_dir, sample_names, [annotator_name])
            continue

        r1_full, r2_full = generate_reads_metagenome(
            genome_files=[available_genomes[taxid] for taxid in genome_ids],
            output_dir=pool_dir,
            amount=[total_by_taxid[taxid] for taxid in genome_ids],
            read_length=read_length,
            sample_name="full",
            annotator_name=annotator_name,
            model=model,
            genome_ids=genome_ids,
            cpus=cpus,
            seed=seed,
        )
        split_metagenome_to_samples(
            r1_full,
            r2_full,
            sample_source_counts,
            output_dir,
            annotator_name,
        )

    shutil.rmtree(pool_dir, ignore_errors=True)


def generate_iss_test_samples(
    genome_dir: str,
    host_genome: str,
    output_dir: str,
    n_samples: int = 10,
    total_reads: int = 2000,
    host_fraction: Union[str, float] = "RANDOM",
    seed: int = 42,
    model: str = "hiseq",
    genomes: Optional[Sequence[str]] = None,
    cpus: int = 2,
) -> List[str]:
    """
    Generate one host pool and one metagenome pool, then split them into samples.
    """
    os.makedirs(output_dir, exist_ok=True)
    rng = random.Random(seed)
    samples = [str(i) for i in range(1, int(n_samples) + 1)]

    if str(host_fraction).upper() == "RANDOM":
        host_counts = {sample: rng.randint(0, int(total_reads)) for sample in samples}
    else:
        host_n = int(float(host_fraction) * int(total_reads))
        host_counts = {sample: host_n for sample in samples}
    meta_counts = {sample: max(0, int(total_reads) - host_counts[sample]) for sample in samples}

    if genomes:
        genome_files = [
            os.path.join(genome_dir, genome) if not os.path.isabs(genome) else genome
            for genome in genomes
        ]
        genome_files = [path for path in genome_files if os.path.exists(path)]
    else:
        genome_files = sorted(
            glob.glob(os.path.join(genome_dir, "*.fa"))
            + glob.glob(os.path.join(genome_dir, "*.fna"))
            + glob.glob(os.path.join(genome_dir, "*.fasta"))
        )

    pool_dir = os.path.join(output_dir, ".iss_full")
    os.makedirs(pool_dir, exist_ok=True)

    host_r1 = os.path.join(pool_dir, "host_R1.fastq")
    host_r2 = os.path.join(pool_dir, "host_R2.fastq")
    total_host = sum(host_counts.values())
    if total_host > 0 and host_genome and os.path.exists(host_genome):
        generate_reads_genome(
            host_genome,
            os.path.join(pool_dir, "host"),
            total_host,
            model=model,
            cpus=cpus,
            seed=seed,
        )
    else:
        _write_empty_fastq_pair(host_r1, host_r2)

    total_meta = sum(meta_counts.values())
    meta_r1 = os.path.join(pool_dir, "meta_full_R1.fastq")
    meta_r2 = os.path.join(pool_dir, "meta_full_R2.fastq")
    if genome_files and total_meta > 0:
        generate_reads_metagenome(
            genome_files=genome_files,
            output_dir=pool_dir,
            amount=[1] * len(genome_files),
            total_amount=total_meta,
            sample_name="meta",
            annotator_name="full",
            model=model,
            cpus=cpus,
            seed=None if seed is None else seed + 1,
        )
    else:
        _write_empty_fastq_pair(meta_r1, meta_r2)

    tmp_host = {
        sample: (
            os.path.join(pool_dir, f"{sample}_host_R1.fastq"),
            os.path.join(pool_dir, f"{sample}_host_R2.fastq"),
        )
        for sample in samples
    }
    tmp_meta = {
        sample: (
            os.path.join(pool_dir, f"{sample}_meta_R1.fastq"),
            os.path.join(pool_dir, f"{sample}_meta_R2.fastq"),
        )
        for sample in samples
    }
    split_paired_fastq_by_counts(host_r1, host_r2, host_counts, tmp_host)
    split_paired_fastq_by_counts(meta_r1, meta_r2, meta_counts, tmp_meta)

    outputs = []
    for sample in samples:
        for mate, host_path, meta_path in (
            ("R1", tmp_host[sample][0], tmp_meta[sample][0]),
            ("R2", tmp_host[sample][1], tmp_meta[sample][1]),
        ):
            dest = os.path.join(output_dir, f"{sample}_full_{mate}.fastq")
            with open(dest, "w") as out:
                wrote = False
                for src in (host_path, meta_path):
                    if not os.path.exists(src):
                        continue
                    with open(src) as handle:
                        text = handle.read()
                    if text.strip():
                        out.write(text)
                        if not text.endswith("\n"):
                            out.write("\n")
                        wrote = True
                if not wrote:
                    out.write("\n")
            outputs.append(dest)

    shutil.rmtree(pool_dir, ignore_errors=True)
    return outputs

def _resolve_r_executable() -> Tuple[str, Optional[str]]:
    """Resolve R binary and optional library path from user/repo config or PATH."""
    cfg = load_config()
    configured = (cfg.get("r_path") or "R").strip() or "R"
    lib = (cfg.get("r_lib_path") or "").strip()
    r_lib_path = lib or None
    resolved = resolve_executable(configured, tool_key="R")
    token = (resolved or configured).split()[0]
    if os.path.isfile(token) and os.access(token, os.X_OK):
        return token, r_lib_path
    found = shutil.which(configured) or shutil.which("R") or shutil.which("Rscript")
    if found:
        return found, r_lib_path
    raise FileNotFoundError(
        f"R executable not found ({configured!r}). "
        "The R generative package is optional; use the Python regenerator "
        "or set SAMOVAR_USE_R=1 after installing R / samovaR."
    )


def _samovar_annotation_regenerate_python(
    annotation_dir: str,
    config_samovar_dict: dict,
    output_dir: str,
) -> None:
    from samovar.regenerate import regenerate_annotation_tables, write_samovar_config_defaults

    cfg = write_samovar_config_defaults(config_samovar_dict)
    regenerate_annotation_tables(annotation_dir, output_dir, cfg)


def samovar_annotation_regenerate(
    annotation_dir: str,
    config_samovar: str = None,
    output_dir: str = None
) -> None:
    """Regenerate taxonomy tables to abundance CSVs.

    Modes (``regeneration_mode`` in config):

    - ``preserve`` (default): observed counts, no generative remodelling.
    - ``glm``: R samovaR boil when ``SAMOVAR_USE_R=1``, else Python glm analog.
    - ``bootstrap``: column bootstrap of observed profiles.
    - ``vae``: latent-factor generative sampling.
    """
    from samovar.regenerate import normalize_regeneration_mode

    tmp_config_path = None
    if config_samovar is None:
        config_samovar_dict = {
            "threshold_amount": 1e-5,
            "plot_log": False,
            "N_reads": 1000,
            "regeneration_mode": "preserve",
        }
    else:
        with open(config_samovar, "r") as f:
            config_samovar_dict = yaml.safe_load(f) or {}

    if output_dir is None:
        output_dir = config_samovar_dict.get("output_dir")
        if not output_dir:
            raise ValueError("output_dir is required")

    mode = normalize_regeneration_mode(config_samovar_dict.get("regeneration_mode", "preserve"))
    use_r = os.environ.get("SAMOVAR_USE_R", "0").strip().lower() in {"1", "true", "yes"}
    if mode != "glm" or not use_r:
        _samovar_annotation_regenerate_python(
            annotation_dir, config_samovar_dict, output_dir
        )
        return

    annotation_regenerate = str(annotation_regenerate_r())
    if not os.path.isfile(annotation_regenerate):
        raise FileNotFoundError(
            f"annotation_regenerate.R not found at {annotation_regenerate}"
        )

    r_path, r_lib_path = _resolve_r_executable()
    os.makedirs(output_dir, exist_ok=True)
    env = os.environ.copy()
    if r_lib_path:
        env["R_LIBS"] = r_lib_path
        env["R_LIBS_USER"] = r_lib_path

    config_path = config_samovar
    if config_path is None:
        handle = tempfile.NamedTemporaryFile(
            mode="w", suffix=".yaml", delete=False, encoding="utf-8"
        )
        yaml.dump(config_samovar_dict, handle)
        handle.close()
        tmp_config_path = handle.name
        config_path = tmp_config_path

    cmd = [
        r_path,
        "--vanilla",
        "-s",
        "-f",
        annotation_regenerate,
        "--args",
        "--config",
        str(config_path),
        "--annotation_dir",
        str(annotation_dir),
        "--output_dir",
        str(output_dir),
    ]
    try:
        subprocess.run(cmd, check=True, env=env)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            "samovar_annotation_regenerate failed while running R "
            f"(exit {exc.returncode}). Command: {' '.join(cmd)}"
        ) from exc
    finally:
        if tmp_config_path:
            try:
                os.unlink(tmp_config_path)
            except OSError:
                pass
