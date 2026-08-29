"""
Module for converting tables to simulated reads using ISS-like functionality.
"""

import os
import re
import glob
import random
import pandas as pd
import subprocess
from contextlib import contextmanager
from .genome_fetcher import fetch_genome, default_entrez_email
from typing import Any, Dict, Iterable, Iterator, List, Optional, Sequence, Tuple, Union
import yaml
import json
import tempfile
import warnings
import shutil
from pathlib import Path
from samovar.paths import annotation_regenerate_r, iss_executable, load_config, user_config_dir
from samovar.seqio import (
    concat_fastq_files,
    fastq_pair_paths,
    find_genome_file,
    gzip_file,
    gunzip_file,
    is_gzip_file,
    is_gzip_path,
    iter_fastq_records,
    list_fasta_files,
    open_text,
    sequence_stem,
    uncompressed_fasta_for_tool,
    write_text_lines,
)


# Thin R CLI used with the optional samovaR package (GitHub branch r-package).
# Not an R package source; it only calls exported samovaR functions.
R_REGENERATE_DRIVER = r'''library(tidyverse)
library(samovaR)

args <- commandArgs(trailingOnly = TRUE)
annotation_dir <- NULL
output_dir <- NULL
config_samovar <- NULL
i <- 1
while (i <= length(args)) {
  if (args[i] == "--annotation_dir") {
    annotation_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] %in% c("--abundance_dir", "--abundance")) {
    annotation_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--output_dir") {
    output_dir <- args[i + 1]; i <- i + 2
  } else if (args[i] == "--config") {
    config_samovar <- args[i + 1]; i <- i + 2
  } else {
    i <- i + 1
  }
}
if (is.null(config_samovar)) {
  config <- list()
} else {
  config <- unpack_config(config_samovar)
}
if (!("N" %in% names(config))) config$N <- 1
if (!("N_reads" %in% names(config))) config$N_reads <- 100
if (!("plot_log" %in% names(config))) config$plot_log <- FALSE
if ("output_dir" %in% names(config)) output_dir <- config$output_dir
if (is.null(annotation_dir)) stop("--annotation_dir is required")
if (!is.null(output_dir)) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
} else {
  output_dir <- "."
}

is_abundance_csv <- function(path) {
  hdr <- tryCatch(colnames(utils::read.csv(path, nrows = 1, check.names = FALSE)), error = function(e) character())
  any(tolower(hdr) %in% c("taxid", "otu", "otu_id", "feature")) &&
    (any(startsWith(hdr, "N_")) || sum(!tolower(hdr) %in% c("taxid", "otu", "otu_id", "feature")) >= 1)
}

load_abundance_list <- function(root) {
  files <- list.files(root, pattern = "\\.csv$", full.names = TRUE)
  files <- files[!grepl("^combined_annotation_table|^\\.", basename(files))]
  out <- list()
  for (f in files) {
    df <- tryCatch(utils::read.csv(f, check.names = FALSE), error = function(e) NULL)
    if (is.null(df) || !nrow(df)) next
    id_idx <- which(tolower(names(df)) %in% c("taxid", "otu", "otu_id", "feature"))
    if (!length(id_idx)) next
    ids <- make.unique(as.character(df[[id_idx[[1]]]]))
    rest <- df[-id_idx[[1]]]
    keep <- vapply(rest, function(x) is.numeric(x) || is.integer(x), logical(1))
    if (!any(keep)) next
    rest <- rest[, keep, drop = FALSE]
    names(rest) <- sub("^N_", "", names(rest))
    mat <- data.matrix(rest)
    rownames(mat) <- ids
    name <- tools::file_path_sans_ext(basename(f))
    out[[name]] <- table2samovar(mat)
  }
  out
}

samovar_data_list <- load_abundance_list(annotation_dir)
if (!length(samovar_data_list)) {
  samovar_data_long <- read_annotation_dir(annotation_dir)
  samovar_data_list <- annotation2samovar(samovar_data_long)
}
for (i in seq_along(samovar_data_list)) {
  annotator <- names(samovar_data_list)[i]
  tryCatch({
    config$samovar_data <- samovar_data_list[[i]]$copy()
    samovar <- do.call(samovar_preprocess, config)
    new_data <- samovar_boil(samovar, N = config$N)
    result_df <- as.data.frame(round(new_data$data * config$N_reads))
    result_df <- tibble::rownames_to_column(result_df, "taxid")
    colnames(result_df)[-1] <- paste0("N_", colnames(result_df)[-1])
    write.csv(
      result_df,
      paste0(output_dir, "/", annotator, ".csv"),
      row.names = FALSE
    )
  }, error = function(e) {
    warning(paste0("samovaR generation failed for ", annotator, ": ", conditionMessage(e)))
  })
}
'''


SKIP_TAXIDS = {"0", "nan", "None", ""}
CONTIG_SPACER = "N" * 500


_ISS_EXTRA_FLAGS: List[str] = []


def _normalize_iss_extra(extra: Optional[Sequence[str]]) -> List[str]:
    from samovar.table_regenerators import extra_flags_argv

    if extra is None:
        return []
    if isinstance(extra, str):
        return extra_flags_argv(extra)
    items = [str(x) for x in extra]
    if items and all(len(x) <= 1 for x in items):
        return extra_flags_argv("".join(items))
    out: List[str] = []
    for item in items:
        if " " in item.strip():
            out.extend(extra_flags_argv(item))
        elif item.strip():
            out.append(item)
    return out


@contextmanager
def iss_cli_extra_flags(extra: Optional[Sequence[str]] = None) -> Iterator[None]:
    """Append import/prepare extra flags to every ISS CLI invocation."""
    global _ISS_EXTRA_FLAGS
    prev = list(_ISS_EXTRA_FLAGS)
    _ISS_EXTRA_FLAGS = _normalize_iss_extra(extra)
    try:
        yield
    finally:
        _ISS_EXTRA_FLAGS = prev


def _iss_cmd(*args: str) -> List[str]:
    return [iss_executable(), *args, *_ISS_EXTRA_FLAGS]


def _run_iss(cmd: Sequence[str]) -> None:
    result = subprocess.run(list(cmd), check=False)
    if result.returncode != 0:
        raise RuntimeError(
            f"ISS command failed (exit {result.returncode}): {' '.join(cmd)}"
        )


def parse_annotation_table(table_path: str) -> pd.DataFrame:
    """Parse a long annotation CSV or an OTU/abundance CSV into ``taxid`` + ``N_*``.

    Long tables (``taxID_*``) are counted per annotator. Phyloseq-style OTU
    tables are rewritten to ``N_<sample>``.
    """
    from samovar.abundance import (
        annotator_name,
        input_to_abundance_tables,
        load_table_input,
        normalize_abundance_table,
        n_sample_columns,
    )

    loaded = load_table_input(table_path)
    tables = input_to_abundance_tables(loaded)
    if not tables:
        return pd.DataFrame(columns=["taxid"])
    if len(tables) == 1:
        return normalize_abundance_table(next(iter(tables.values())))
    result = pd.DataFrame()
    for name, table in tables.items():
        n_cols = n_sample_columns(table)
        if not n_cols:
            continue
        counts = table[["taxid"]].copy()
        counts[f"N_{annotator_name(name)}"] = table[n_cols].sum(axis=1)
        if result.empty:
            result = counts
        else:
            result = result.merge(counts, on="taxid", how="outer")
    if result.empty:
        return pd.DataFrame(columns=["taxid"])
    result = result.fillna(0)
    result["taxid"] = result["taxid"].astype(str)
    return result

def get_genome_file(genome_dir: str, taxid: str) -> str:
    """
    Get the path to a genome file, checking multiple possible extensions
    (processed, gzipped, and uncompressed FASTA).
    """
    return find_genome_file(genome_dir, taxid)


def _n_columns(df: pd.DataFrame) -> List[str]:
    from samovar.abundance import n_sample_columns

    cols = n_sample_columns(df)
    if cols:
        return cols
    return [col for col in df.columns if col != "taxid" and str(col).startswith("N_")]


def _annotator_from_n_col(column: str) -> str:
    match = re.search(r"N_(.*?)(?:_[0-9]*)?$", column)
    return match.group(1) if match else "any"


def _safe_fasta_id(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9.-]", "_", str(value))
    return cleaned or "genome"


def _genome_id_from_path(path: str) -> str:
    return _safe_fasta_id(sequence_stem(path))


def _open_text(path: str):
    return open_text(path)


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
        target = max(0, int(total_amount))
        raw = [n / total * target for n in scaled]
        scaled = [int(value) for value in raw]
        remainder = target - sum(scaled)
        order = sorted(
            range(len(raw)),
            key=lambda i: (raw[i] - scaled[i], -i),
            reverse=True,
        )
        for i in order[:remainder]:
            scaled[i] += 1
    return scaled


def _iss_reads_for_pairs(n_pairs: int) -> int:
    """ISS counts individual reads; SamovaR abundance rows count read pairs."""
    return max(0, int(n_pairs)) * 2


def _write_empty_fastq_pair(r1: str, r2: str) -> None:
    for path in (r1, r2):
        os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
        with open_text(path, "wt"):
            pass


def _iss_fastq_candidates(prefix: str, mate: str) -> List[str]:
    return [
        f"{prefix}_{mate}.fastq",
        f"{prefix}_{mate}.fastq.gz",
        f"{prefix}_{mate}.fq",
        f"{prefix}_{mate}.fq.gz",
        f"{prefix}.{mate}.fastq",
        f"{prefix}.{mate}.fastq.gz",
    ]


def _first_existing_fastq(paths: Sequence[str]) -> Optional[str]:
    nonempty = [p for p in paths if os.path.exists(p) and os.path.getsize(p) > 0]
    if nonempty:
        return nonempty[0]
    empty = [p for p in paths if os.path.exists(p)]
    return empty[0] if empty else None


def _fastq_has_records(path: str) -> bool:
    from samovar.seqio import fastq_has_reads

    return fastq_has_reads(path)


def ensure_uncompressed_iss_pair(prefix: str) -> Tuple[str, str]:
    """Copy/gunzip ISS output to ``{prefix}_R1.fastq`` / ``_R2.fastq``."""
    outs: List[str] = []
    for mate in ("R1", "R2"):
        dest = f"{prefix}_{mate}.fastq"
        found = _first_existing_fastq(_iss_fastq_candidates(prefix, mate))
        if found is None:
            outs.append(dest)
            continue
        os.makedirs(os.path.dirname(dest) or ".", exist_ok=True)
        if is_gzip_file(found):
            gunzip_file(found, dest=dest, remove_source=is_gzip_path(found))
        elif os.path.abspath(found) != os.path.abspath(dest):
            shutil.copy2(found, dest)
        outs.append(dest)
    return outs[0], outs[1]


def _cleanup_iss_tmp(prefix: str) -> None:
    directory = os.path.dirname(prefix) or "."
    base = os.path.basename(prefix)
    for path in glob.glob(os.path.join(directory, f"{base}*")):
        name = os.path.basename(path)
        if name.endswith("_R1.fastq") or name.endswith("_R2.fastq"):
            continue
        if name.endswith("_R1.fastq.gz") or name.endswith("_R2.fastq.gz"):
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
    rec1_iter = iter_fastq_records(r1_path)
    rec2_iter = iter_fastq_records(r2_path)
    for rec1, rec2 in zip(rec1_iter, rec2_iter):
        if not rec1[0].strip():
            continue
        yield list(rec1), list(rec2)


def _write_fastq_records(path: str, records: Iterable[List[str]]) -> None:
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    write_text_lines(path, ("".join(record) for record in records))


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
    if not os.path.exists(r1_path):
        for sample, (r1_out, r2_out) in output_paths.items():
            _write_empty_fastq_pair(r1_out, r2_out)
        return {sample: 0 for sample in output_paths}

    if not os.path.exists(r2_path) or not _fastq_has_records(r2_path):
        recs = [list(rec) for rec in iter_fastq_records(r1_path) if rec[0].strip()]
        empty = ["", "", "", ""]
        pairs = [(rec, list(empty)) for rec in recs]
    else:
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
    gzip_reads: bool = False,
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
        r1_out, r2_out = fastq_pair_paths(
            os.path.join(output_dir, f"{sample}_{annotator_name}"),
            gzip_reads=gzip_reads,
        )
        os.makedirs(output_dir, exist_ok=True)
        output_handles[sample] = (open_text(r1_out, "wt"), open_text(r2_out, "wt"))

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
                handle.close()

    return written


def _maybe_gzip_fastq_pair(r1: str, r2: str, gzip_reads: bool) -> Tuple[str, str]:
    if not gzip_reads:
        return r1, r2
    out = []
    for path in (r1, r2):
        if path and os.path.exists(path) and not is_gzip_path(path):
            out.append(str(gzip_file(path)))
        else:
            out.append(path)
    return out[0], out[1]


def generate_reads_genome(
    genome_file: str,
    output_file: str,
    amount: int,
    read_length: int = 150,
    model: str = "hiseq",
    cpus: int = 2,
    seed: Optional[int] = None,
    gzip_reads: bool = False,
) -> Tuple[str, str]:
    """
    Generate ``amount`` paired records from a genome file.

    ISS requires an uncompressed FASTA path; gzipped genomes are decompressed
    to a temp file for the ISS call only.
    """
    r1 = f"{output_file}_R1.fastq"
    r2 = f"{output_file}_R2.fastq"
    if genome_file is not None and os.path.exists(genome_file):
        os.makedirs(os.path.dirname(output_file) or ".", exist_ok=True)
        with uncompressed_fasta_for_tool(genome_file) as plain_fasta:
            cmd = _iss_cmd(
                "generate",
                "--genomes", plain_fasta,
                "--model", model,
                "--output", output_file,
                "--n_reads", str(_iss_reads_for_pairs(amount)),
                "--cpus", str(cpus),
            )
            if seed is not None:
                cmd.extend(["--seed", str(seed)])
            _run_iss(cmd)
        _cleanup_iss_tmp(output_file)
        ensure_uncompressed_iss_pair(output_file)
        if amount > 0 and not _fastq_has_records(r1):
            raise RuntimeError(f"ISS did not produce reads for {output_file}")
    return _maybe_gzip_fastq_pair(r1, r2, gzip_reads)


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
    gzip_reads: bool = False,
) -> Tuple[str, str]:
    """
    Generate paired records from a metagenome in a single ISS invocation.

    ``amount`` and ``total_amount`` count pairs. ISS's readcount file counts
    individual reads, so each value is doubled at the CLI boundary. Combined
    genomes are written uncompressed for ISS, then gzipped afterwards.
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
        return _maybe_gzip_fastq_pair(output_file_R1, output_file_R2, gzip_reads)

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
        [_iss_reads_for_pairs(amount_by_id[gid]) for gid in written_ids],
    )
    if (
        total_reads <= 0
        or not os.path.exists(combined_fasta)
        or os.path.getsize(combined_fasta) == 0
    ):
        _write_empty_fastq_pair(output_file_R1, output_file_R2)
        return _maybe_gzip_fastq_pair(output_file_R1, output_file_R2, gzip_reads)

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
    ensure_uncompressed_iss_pair(output_prefix)
    if os.path.exists(combined_fasta):
        try:
            gzip_file(combined_fasta)
        except OSError:
            pass

    if not _fastq_has_records(output_file_R1) or not os.path.exists(output_file_R2):
        raise RuntimeError(
            f"ISS did not produce reads for {output_prefix}"
        )

    return _maybe_gzip_fastq_pair(output_file_R1, output_file_R2, gzip_reads)

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
    gzip_reads: bool = False,
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
        gzip_reads=gzip_reads,
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
    gzip_genomes: bool = True,
    reannotation_level: str = "taxid",
    max_genome_mb: Any = None,
    genome_skip_list: Any = None,
) -> Dict[str, str]:
    from samovar.genome_resolve import resolve_genome_file

    available = {}
    from samovar.regenerate import finite_max_genomes

    limit = finite_max_genomes(max_genomes, default_from_env=False)
    for taxid in taxids:
        if limit is not None and len(available) >= limit:
            break
        taxid = str(taxid).split(".")[0]
        if taxid in SKIP_TAXIDS or taxid in available:
            continue
        genome_file = resolve_genome_file(
            taxid,
            genome_dir,
            email,
            level=reannotation_level,
            reference_only=reference_only,
            gzip_genomes=gzip_genomes,
            max_genome_mb=max_genome_mb,
            genome_skip_list=genome_skip_list,
        )
        if genome_file is not None:
            available[taxid] = genome_file
    return available


def _emit_empty_for_annotators(
    output_dir: str,
    sample_names: Sequence[str],
    annotator_names: Sequence[str],
    gzip_reads: bool = False,
) -> None:
    os.makedirs(output_dir, exist_ok=True)
    for sample_name in sample_names:
        for annotator_name in annotator_names:
            r1, r2 = fastq_pair_paths(
                os.path.join(output_dir, f"{sample_name}_{annotator_name}"),
                gzip_reads=gzip_reads,
            )
            _write_empty_fastq_pair(r1, r2)


def process_abundance_table(
    table: Union[str, pd.DataFrame],
    genome_dir: str,
    output_dir: str,
    total_amount: int = None,
    email: Optional[str] = None,
    reference_only: bool = True,
    model: str = "hiseq",
    read_length: int = 150,
    sample_name: str = None,
    gzip_genomes: bool = True,
    gzip_reads: bool = False,
    reannotation_level: str = "taxid",
    max_genome_mb: Any = None,
    genome_skip_list: Any = None,
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
    from samovar.abundance import normalize_abundance_table

    abundance_table = normalize_abundance_table(abundance_table)
    
    if sample_name is None:
        if isinstance(table, str):
            sample_name = os.path.basename(table).split(".")[0]
        else:
            sample_name = "merged"

    email = email or default_entrez_email()

    if "taxid" not in abundance_table.columns:
        _emit_empty_for_annotators(output_dir, [sample_name], ["any"], gzip_reads=gzip_reads)
        warnings.warn("Abundance table has no taxid column; emitted empty FASTQ files")
        return None

    abundance_table = abundance_table.copy()
    abundance_table["taxid"] = abundance_table["taxid"].astype(str).str.split(".").str[0]
    available_genomes = _resolve_genomes_for_taxids(
        abundance_table["taxid"],
        genome_dir,
        email,
        reference_only,
        gzip_genomes=gzip_genomes,
        reannotation_level=reannotation_level,
        max_genome_mb=max_genome_mb,
        genome_skip_list=genome_skip_list,
    )

    N_cols = _n_columns(abundance_table)
    annotators = [_annotator_from_n_col(col) for col in N_cols] or ["any"]

    # Gracefully handle all-unclassified/unknown samples.
    if not available_genomes:
        _emit_empty_for_annotators(output_dir, [sample_name], annotators, gzip_reads=gzip_reads)
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
            gzip_reads=gzip_reads,
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
    abundance_dir: Optional[str] = None,
    regeneration_config: Optional[dict] = None,
    gzip_genomes: bool = True,
    gzip_reads: bool = False,
    reannotation_level: str = "taxid",
    max_genome_mb: Any = None,
    genome_skip_list: Any = None,
) -> None:
    """Simulate FASTQ from abundance tables (abundance2iss).

    ``annotation2iss`` is annotation2abundance then this function. Long
    ``*.annotation.csv`` / OTU tables are converted to ``taxid`` + ``N_*``
    CSVs first; generative table regen belongs to the ``regenerate_tables``
    checkpoint (or a non-direct mode in ``regeneration_config`` when dest is empty).
    """
    from samovar.abundance import (
        convert_to_abundance_dir,
        has_abundance_tables,
    )
    from samovar.regenerate import (
        coerce_seed,
        write_samovar_config_defaults,
    )

    regen_cfg = write_samovar_config_defaults(dict(regeneration_config or {}))
    from samovar.regenerate import finite_max_genomes

    if max_genomes is None:
        max_genomes = regen_cfg.get("max_genomes")
    max_genomes = finite_max_genomes(max_genomes, default_from_env=False)
    if max_genome_mb is None:
        max_genome_mb = regen_cfg.get("max_genome_mb")
    if genome_skip_list is None:
        genome_skip_list = regen_cfg.get("genome_skip_list")
    if seed is None:
        seed = regen_cfg.get("seed")
    seed = coerce_seed(seed)

    table_paths = list(table_paths)
    email = email or default_entrez_email()
    if sample_names is None:
        sample_names = [os.path.basename(path).split(".")[0] for path in table_paths]
    else:
        sample_names = [
            name or os.path.basename(path).split(".")[0]
            for name, path in zip(sample_names, table_paths)
        ]

    if not annotation_dir and table_paths:
        parent = os.path.dirname(os.path.abspath(table_paths[0]))
        if all(os.path.dirname(os.path.abspath(p)) == parent for p in table_paths) and any(
            os.path.basename(p).endswith(".annotation.csv") for p in table_paths
        ):
            annotation_dir = parent

    sample_tables: Dict[str, pd.DataFrame] = {}
    regen_dir = os.path.join(output_dir, ".regenerated_abundance")
    abund = str(
        abundance_dir
        or regen_cfg.get("abundance_dir")
        or regen_dir
    )
    observed = regen_cfg.get("observed_abundance_dir") or ""
    source = ""
    if observed and has_abundance_tables(observed):
        source = observed
    elif annotation_dir:
        source = annotation_dir
    elif regen_cfg.get("annotation_dir"):
        source = str(regen_cfg.get("annotation_dir"))
    elif table_paths:
        source = os.path.dirname(os.path.abspath(table_paths[0]))
    if not has_abundance_tables(abund) and source:
        convert_to_abundance_dir(source, abund, regen_cfg)
    existing = _sample_tables_from_abundance_dir(abund, sample_names)
    if existing:
        sample_tables = existing
        sample_names = list(sample_tables.keys())
    elif table_paths:
        for path, sample_name in zip(table_paths, sample_names):
            table = parse_annotation_table(path)
            if not table.empty and "taxid" in table.columns:
                table = table.copy()
                table["taxid"] = table["taxid"].astype(str).str.split(".").str[0]
            sample_tables[sample_name] = table

    if not sample_tables:
        _emit_empty_for_annotators(
            output_dir, sample_names or ["1"], ["any"], gzip_reads=gzip_reads
        )
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
        gzip_genomes=gzip_genomes,
        reannotation_level=reannotation_level,
        max_genome_mb=max_genome_mb,
        genome_skip_list=genome_skip_list,
    )
    if not available_genomes:
        _emit_empty_for_annotators(output_dir, sample_names, annotators, gzip_reads=gzip_reads)
        warnings.warn("No genome files available for any taxid; emitted empty FASTQ files")
        return

    from samovar.reads_generators import simulate_from_sample_tables

    annotator_map: Dict[str, str] = dict(annotator_cols)
    simulate_from_sample_tables(
        sample_tables,
        available_genomes,
        output_dir,
        annotator_cols=annotator_map,
        total_amount=total_amount,
        model=model,
        read_length=read_length,
        cpus=cpus,
        seed=int(seed or 42),
        gzip_reads=gzip_reads,
        reads_generator=regen_cfg.get("reads_generator") or regen_cfg.get("reads-generator"),
        extra_flags=None,
        config=regen_cfg,
        genome_dir=genome_dir,
    )


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
    extra_flags: Optional[Sequence[str]] = None,
    abundance_table: Optional[Union[str, pd.DataFrame]] = None,
    max_genomes: Any = None,
) -> List[str]:
    """
    Generate one host pool and one metagenome pool, then split them into samples.
    """
    from samovar.regenerate import coerce_seed
    from samovar.seqio import taxid_from_fasta_name

    os.makedirs(output_dir, exist_ok=True)
    seed = coerce_seed(seed)
    has_table = isinstance(abundance_table, pd.DataFrame) or (
        abundance_table not in (None, "", False)
    )
    if has_table:
        with iss_cli_extra_flags(extra_flags):
            return generate_iss_from_abundance_table(
                abundance_table,
                genome_dir=genome_dir,
                output_dir=output_dir,
                host_genome=host_genome,
                host_fraction=host_fraction,
                total_reads=total_reads,
                seed=seed,
                model=model,
                genomes=genomes,
                cpus=cpus,
                max_genomes=max_genomes,
            )
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
        genome_files = [str(p) for p in list_fasta_files(genome_dir, nucleotide=True, protein=False)]
    from samovar.regenerate import cap_sequence

    genome_files = cap_sequence(genome_files, max_genomes)

    pool_dir = os.path.join(output_dir, ".iss_full")
    os.makedirs(pool_dir, exist_ok=True)

    with iss_cli_extra_flags(extra_flags):
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
            host_tax = taxid_from_fasta_name(host_genome) or ""
            if host_tax:
                from samovar.camisim import tag_fastq_file

                tag_fastq_file(host_r1, host_r1, "", str(host_tax))
                if os.path.exists(host_r2):
                    tag_fastq_file(host_r2, host_r2, "", str(host_tax))
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
            concat_fastq_files([host_path, meta_path], dest)
            outputs.append(dest)

    if int(total_reads) > 0 and not any(_fastq_has_records(path) for path in outputs):
        raise RuntimeError(
            f"ISS generate wrote no FASTQ records under {output_dir} "
            f"(n_samples={n_samples}, total_reads={total_reads})"
        )

    shutil.rmtree(pool_dir, ignore_errors=True)
    return outputs


def generate_iss_from_abundance_table(
    table: Union[str, pd.DataFrame],
    genome_dir: str,
    output_dir: str,
    host_genome: str = "",
    host_fraction: Union[str, float] = 0,
    total_reads: Optional[int] = None,
    seed: int = 42,
    model: str = "hiseq",
    genomes: Optional[Sequence[str]] = None,
    cpus: int = 2,
    annotator_name: str = "full",
    max_genomes: Any = None,
) -> List[str]:
    """ISS from a taxid × ``N_<sample>`` table; names ``{sample}_{annotator}_R*``."""
    from samovar.seqio import taxid_from_fasta_name

    annotator_name = str(annotator_name or "full")

    if isinstance(table, str):
        abundance = pd.read_csv(table)
    else:
        abundance = table.copy()
    from samovar.abundance import normalize_abundance_table

    abundance = normalize_abundance_table(abundance)
    from samovar.regenerate import cap_abundance_table, finite_max_genomes

    abundance = cap_abundance_table(abundance, max_genomes)
    if "taxid" not in abundance.columns:
        raise ValueError("abundance table must have a taxid column")
    abundance["taxid"] = abundance["taxid"].astype(str).str.split(".").str[0]
    n_cols = _n_columns(abundance)
    sample_names = []
    for col in n_cols:
        name = str(col)
        if name.lower().startswith("n_"):
            name = name[2:]
        sample_names.append(name or str(col))
    if not n_cols:
        _write_empty_fastq_pair(
            os.path.join(output_dir, f"1_{annotator_name}_R1.fastq"),
            os.path.join(output_dir, f"1_{annotator_name}_R2.fastq"),
        )
        return [
            os.path.join(output_dir, f"1_{annotator_name}_R1.fastq"),
            os.path.join(output_dir, f"1_{annotator_name}_R2.fastq"),
        ]

    sample_tables: Dict[str, pd.DataFrame] = {}
    for sample, col in zip(sample_names, n_cols):
        piece = abundance[["taxid", col]].copy()
        piece = piece.rename(columns={col: f"N_{sample}"})
        sample_tables[str(sample)] = piece

    if genomes:
        from samovar.regenerate import cap_sequence

        genomes = cap_sequence(list(genomes), max_genomes)
        genome_files = {
            _genome_id_from_path(p): (
                os.path.join(genome_dir, p) if not os.path.isabs(p) else p
            )
            for p in genomes
        }
        available = {k: v for k, v in genome_files.items() if os.path.exists(v)}
    else:
        available = _resolve_genomes_for_taxids(
            abundance["taxid"].astype(str).tolist(),
            genome_dir,
            default_entrez_email(),
            True,
            max_genomes=finite_max_genomes(max_genomes),
        )
    outputs: List[str] = []
    pool_dir = os.path.join(output_dir, ".iss_full")
    os.makedirs(pool_dir, exist_ok=True)
    sample_source_counts: Dict[str, Dict[str, int]] = {}
    total_by_taxid: Dict[str, int] = {}
    for sample_name, tbl in sample_tables.items():
        counts: Dict[str, int] = {}
        col = _n_columns(tbl)[0]
        amounts = tbl[col].tolist()
        if total_reads:
            amounts = _scale_amounts(amounts, int(total_reads))
        for taxid, n_reads in zip(tbl["taxid"].astype(str), amounts):
            if taxid not in available or int(n_reads) <= 0:
                continue
            counts[taxid] = counts.get(taxid, 0) + int(n_reads)
            total_by_taxid[taxid] = total_by_taxid.get(taxid, 0) + int(n_reads)
        sample_source_counts[sample_name] = counts
    genome_ids = [t for t, n in total_by_taxid.items() if n > 0]
    if genome_ids:
        r1_full, r2_full = generate_reads_metagenome(
            genome_files=[available[t] for t in genome_ids],
            output_dir=pool_dir,
            amount=[total_by_taxid[t] for t in genome_ids],
            sample_name="pool",
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
            gzip_reads=False,
        )
        # split writes {sample}_full_R* already because annotator_name is full.
    else:
        for sample_name in sample_tables:
            _write_empty_fastq_pair(
                os.path.join(output_dir, f"{sample_name}_{annotator_name}_R1.fastq"),
                os.path.join(output_dir, f"{sample_name}_{annotator_name}_R2.fastq"),
            )
    if host_genome and os.path.exists(host_genome) and str(host_fraction).upper() not in {"0", "0.0", ""}:
        # Optional host spike-in on top of the table (same fractions as generate).
        host_tax = taxid_from_fasta_name(host_genome) or ""
        rng = random.Random(seed)
        for sample_name in sample_tables:
            dest_r1 = os.path.join(output_dir, f"{sample_name}_{annotator_name}_R1.fastq")
            dest_r2 = os.path.join(output_dir, f"{sample_name}_{annotator_name}_R2.fastq")
            if str(host_fraction).upper() == "RANDOM":
                host_n = rng.randint(0, int(total_reads or 0))
            else:
                host_n = int(float(host_fraction) * int(total_reads or 0))
            if host_n <= 0:
                continue
            tmp_prefix = os.path.join(pool_dir, f"{sample_name}_host")
            generate_reads_genome(
                host_genome, tmp_prefix, host_n, model=model, cpus=cpus, seed=seed
            )
            hr1, hr2 = f"{tmp_prefix}_R1.fastq", f"{tmp_prefix}_R2.fastq"
            if host_tax:
                from samovar.camisim import tag_fastq_file

                tag_fastq_file(hr1, hr1, "", str(host_tax))
                if os.path.exists(hr2):
                    tag_fastq_file(hr2, hr2, "", str(host_tax))
            concat_fastq_files([hr1, dest_r1], dest_r1 + ".mix")
            shutil.move(dest_r1 + ".mix", dest_r1)
            if os.path.exists(hr2) and os.path.exists(dest_r2):
                concat_fastq_files([hr2, dest_r2], dest_r2 + ".mix")
                shutil.move(dest_r2 + ".mix", dest_r2)
    shutil.rmtree(pool_dir, ignore_errors=True)
    for sample_name in sample_tables:
        outputs.extend(
            [
                os.path.join(output_dir, f"{sample_name}_{annotator_name}_R1.fastq"),
                os.path.join(output_dir, f"{sample_name}_{annotator_name}_R2.fastq"),
            ]
        )
    return outputs


def _resolve_r_executable() -> Tuple[str, Optional[str]]:
    """Resolve R binary and optional library path from user/repo config or PATH."""
    cfg = load_config()
    configured = (cfg.get("r_path") or "R").strip() or "R"
    lib = (cfg.get("r_lib_path") or "").strip()
    r_lib_path = lib or None
    token = configured.split()[0]
    if os.path.isfile(token) and os.access(token, os.X_OK):
        return token, r_lib_path
    found = shutil.which(configured)
    if found:
        return found, r_lib_path
    # Only the default lookup may fall back to Rscript. An explicit invalid
    # r_path is a configuration error and must not be silently ignored.
    if configured == "R":
        found = shutil.which("Rscript")
        if found:
            return found, r_lib_path
    raise FileNotFoundError(
        f"R executable not found ({configured!r}). "
        "Mode 'samovar' needs R plus the optional samovaR regenerator. "
        "Install them separately and set r_path / SAMOVAR_R_REGENERATE, "
        "or use regeneration_mode=direct|bootstrap|vae|glm (Python)."
    )


def _samovar_annotation_regenerate_python(
    annotation_dir: str,
    config_samovar_dict: dict,
    output_dir: str,
) -> None:
    from samovar.regenerate import regenerate_annotation_tables, write_samovar_config_defaults

    cfg = write_samovar_config_defaults(config_samovar_dict)
    regenerate_annotation_tables(annotation_dir, output_dir, cfg)


def _optional_r_regenerator_script() -> Path:
    script = annotation_regenerate_r()
    if script is not None:
        path = Path(script)
        if path.is_file():
            return path
        raise FileNotFoundError(
            "regeneration_mode='samovar' requires the optional R regenerator. "
            "Install samovaR from https://github.com/ctlab/samovar/tree/r-package "
            "(./install.sh R-package) and set SAMOVAR_R_REGENERATE or "
            "annotation_regenerate_r. "
            f"Looked up: {path}"
        )
    dest = user_config_dir() / "annotation_regenerate.R"
    dest.parent.mkdir(parents=True, exist_ok=True)
    if not dest.is_file():
        dest.write_text(R_REGENERATE_DRIVER, encoding="utf-8")
    return dest


def samovar_annotation_regenerate(
    annotation_dir: str,
    config_samovar: str = None,
    output_dir: str = None
) -> None:
    """Regenerate taxonomy tables to abundance CSVs.

    Modes (``regeneration_mode`` / ``table_reads_generator`` in config):

    - ``direct`` (default; alias ``preserve``): observed counts, same samples.
    - ``bootstrap`` / ``vae`` / ``glm``: Python generative models.
    - ``samovar``: optional R regenerator (not part of the Python install).
    - any other name: imported ``table_reads_generator`` (``samovar tools import --type table``).
    """
    from samovar.regenerate import resolve_regeneration_mode

    tmp_config_path = None
    if config_samovar is None:
        config_samovar_dict = {
            "threshold_amount": 1e-5,
            "plot_log": False,
            "N_reads": 1000,
            "regeneration_mode": "direct",
        }
    else:
        with open(config_samovar, "r") as f:
            config_samovar_dict = yaml.safe_load(f) or {}

    if output_dir is None:
        output_dir = config_samovar_dict.get("output_dir")
        if not output_dir:
            raise ValueError("output_dir is required")

    kind, mode = resolve_regeneration_mode(
        config_samovar_dict.get("table_reads_generator")
        or config_samovar_dict.get("regeneration_mode", "direct")
    )
    if not (kind == "builtin" and mode == "samovar"):
        _samovar_annotation_regenerate_python(
            annotation_dir, config_samovar_dict, output_dir
        )
        return

    annotation_regenerate = str(_optional_r_regenerator_script())

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
