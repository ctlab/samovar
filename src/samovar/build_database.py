#!/usr/bin/env python3

import os
import shutil
import subprocess
import logging
import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
from samovar.paths import command_argv

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
from .fasta_processor import process_fasta_directories
from .fasta_processor import preprocess_fasta
from .genome_fetcher import fetch_gff, fetch_proteome, default_entrez_email
from .seqio import (
    gzip_file,
    gunzip_file,
    is_gzip_path,
    iter_seqio_fasta,
    open_text,
    sibling_candidates,
)

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def run_command(cmd: List[str], check: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and return the result.
    
    Args:
        cmd (List[str]): List of command and arguments to execute.
        check (bool, optional): If True, raises CalledProcessError if command fails. Defaults to True.
    
    Returns:
        subprocess.CompletedProcess: Object containing the result of the command execution.
    
    Raises:
        subprocess.CalledProcessError: If check=True and command returns non-zero exit code.
    """
    argv = [str(part) for part in cmd]
    if argv:
        from samovar.paths import resolve_executable

        resolved = command_argv(resolve_executable(argv[0]) or "")
        if resolved and resolved[0]:
            argv[0] = resolved[0]
    logger.info(f"Running command: {' '.join(argv)}")
    return subprocess.run(argv, check=check, text=True)

def get_taxonomy_db(
    db_path: str="kraken_db",
    taxonomy_path: str=None
) -> None:
    """Place NCBI ``nodes.dmp`` / ``names.dmp`` into ``db_path`` as real files.

    Source is ``taxonomy_path`` or the install ``genomes.taxdump`` directory
    (downloaded there once if missing). Copies the dumps so a moved DB folder
    still has ``nodes.dmp``. Replaces leftover symlinks from older builds.
    """
    from samovar.taxdump import ensure_taxdump, link_taxdump_into

    os.makedirs(db_path, exist_ok=True)
    dest = Path(db_path)
    existing_nodes = dest / "nodes.dmp"
    existing_nested = dest / "taxonomy" / "nodes.dmp"

    def _portable(path: Path) -> bool:
        try:
            return path.is_file() and not path.is_symlink()
        except OSError:
            return False

    if _portable(existing_nodes) or _portable(existing_nested):
        logger.info(f"Taxonomy already present under {db_path}; skipping download")
        return

    source = taxonomy_path
    if not source:
        source = str(ensure_taxdump())
    link_taxdump_into(dest, source)

def process_fasta_kraken2(input_file: str, 
                  taxid: str) -> str:
    """Process FASTA file to add proper taxonomy IDs to headers for Kraken2 compatibility.
    
    This function reads a FASTA file and modifies sequence headers to include taxonomy IDs
    in the format required by Kraken2: "seq{index}|kraken:taxid|{taxid}".
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to be added to sequence headers.
    
    Returns:
        str: Path to the temporary processed FASTA file.
    
    Note:
        The function creates a temporary file that should be cleaned up by the caller.
    """
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fa', delete=False)
    records = []
    for i, record in enumerate(iter_seqio_fasta(input_file)):
        record.id = f"seq{i}|kraken:taxid|{taxid}"
        record.description = ""
        records.append(record)
    SeqIO.write(records, temp_fasta.name, "fasta")
    temp_fasta.close()
    return temp_fasta.name

def kraken2_db_workdir(db_path: str) -> Tuple[str, Optional[str]]:
    """Directory Kraken2 scripts may see, and the final path when it differs.

    ``mask_low_complexity.sh`` and ``add_to_library.sh`` expand ``$filename``
    without quotes. A database path that contains whitespace is built in a
    temporary directory and moved into place afterwards. Quoting SAMOVAR's
    own ``--configfile`` does not reach those scripts.
    """
    text = str(db_path)
    if not any(ch.isspace() for ch in text):
        return text, None
    parent = tempfile.gettempdir()
    if any(ch.isspace() for ch in parent) or not os.path.isdir(parent):
        parent = "/tmp"
    stage = tempfile.mkdtemp(prefix="samovar-kraken2-", dir=parent)
    return stage, text


def install_kraken2_db(stage: str, dest: str) -> None:
    """Move a database built without whitespace into ``dest``."""
    destination = Path(dest)
    if destination.exists():
        shutil.rmtree(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.move(stage, destination)


def add_database_kraken2(
    input_file: str,
    taxid: str,
    db_path: str="kraken_db" 
) -> None:
    """Add sequences to the Kraken2 database library.
    
    This function processes a FASTA file and adds its sequences to the Kraken2 database library
    with proper taxonomy IDs. The sequences must be in FASTA format.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to associate with the sequences.
        db_path (str, optional): Path to the Kraken2 database. Defaults to "kraken_db".
    
    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    # Create database directory if it doesn't exist
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    else:
        temp_fasta = process_fasta_kraken2(input_file, taxid)

    try:
        add_cmd = [
            "kraken2-build",
            "--add-to-library",
            temp_fasta,
            "--db",
            db_path
        ]
        run_command(add_cmd)
    finally:
        try:
            os.remove(temp_fasta)
        except OSError:
            pass

    logger.info(f"Database successfully added to {db_path}")


def build_database_kraken2(
    db_path: str="kraken_db", 
    threads: int = 1,
    kmer_len: int = 35,
    minimizer_len: int = 31,
    minimizer_spaces: int = 7,
    skip_maps: bool = True
) -> None:
    """Build a custom Kraken2 database from the library of sequences.
    
    This function builds a Kraken2 database using the sequences that have been added to the library.
    It creates the k-mer database and performs necessary indexing for efficient sequence classification.
    
    Args:
        db_path (str, optional): Path to the database directory. Defaults to "kraken_db".
        threads (int, optional): Number of threads to use for building. Defaults to 1.
        kmer_len (int, optional): Length of k-mers to use. Defaults to 35.
        minimizer_len (int, optional): Length of minimizers to use. Defaults to 31.
        minimizer_spaces (int, optional): Number of spaces in minimizer. Defaults to 7.
        skip_maps (bool, optional): Whether to skip downloading accession maps. Defaults to True.
    
    Note:
        The function will clean up intermediate files after building the database.
        Building a large database can be memory-intensive and time-consuming.
    """
    
    # Build the database
    build_cmd = [
        "kraken2-build",
        "--build",
        "--db", db_path,
        "--threads", str(threads),
        "--kmer-len", str(kmer_len),
        "--minimizer-len", str(minimizer_len),
        "--minimizer-spaces", str(minimizer_spaces),
    ]
    # `--skip-maps` is not available in all kraken2-build versions.
    if skip_maps:
        try:
            help_out = subprocess.run(
                ["kraken2-build", "--help"], check=True, capture_output=True, text=True
            ).stdout
            if "--skip-maps" in help_out:
                build_cmd.append("--skip-maps")
        except Exception:
            pass
    run_command(build_cmd)
    
    # Clean up intermediate files
    clean_cmd = ["kraken2-build", "--clean", "--db", db_path]
    run_command(clean_cmd)
    
    logger.info(f"Database successfully built at {db_path}")


_NUCLEOTIDE_SUFFIXES = (
    ".fna", ".fa", ".fasta", ".fsa", ".ffn",
    ".fna.gz", ".fa.gz", ".fasta.gz", ".fsa.gz", ".ffn.gz",
)


def _is_nucleotide_genome(path: str) -> bool:
    return str(path).lower().endswith(_NUCLEOTIDE_SUFFIXES)


def write_centrifuge_inputs(
    pairs: List[Tuple[str, str]],
    library_fa: str,
    conversion_table: str,
) -> int:
    """Concatenate nucleotide FASTA and write a seqid-to-taxid table.

    ``centrifuge-build`` reads the conversion table as ``seqid<TAB>taxid``.
    Sequence ids are the sanitized FASTA ids, made unique across genomes.
    """
    written = 0
    seen = set()
    with open(library_fa, "w", encoding="utf-8") as fasta, open(
        conversion_table, "w", encoding="utf-8"
    ) as table:
        for path, taxid in pairs:
            tax = str(taxid).split(".")[0]
            for record in iter_seqio_fasta(path):
                base = _sanitize_seqid(record.id) or f"seq{written}"
                seqid = base
                suffix = 2
                while seqid in seen:
                    seqid = f"{base}_{suffix}"
                    suffix += 1
                seen.add(seqid)
                fasta.write(f">{seqid}\n{record.seq}\n")
                table.write(f"{seqid}\t{tax}\n")
                written += 1
    return written


def build_database_centrifuge(db_path: str = "centrifuge_db", threads: int = 1) -> None:
    """Build a Centrifuge index from ``library.fna`` and ``taxonomy/*.dmp``."""
    threads = max(1, int(threads))
    tax_dir = os.path.join(db_path, "taxonomy")
    cmd = [
        "centrifuge-build",
        "-p", str(threads),
        "--conversion-table", os.path.join(db_path, "seqid2taxid.map"),
        "--taxonomy-tree", os.path.join(tax_dir, "nodes.dmp"),
        "--name-table", os.path.join(tax_dir, "names.dmp"),
        os.path.join(db_path, "library.fna"),
        os.path.join(db_path, "centrifuge"),
    ]
    run_command(cmd)
    logger.info(f"Centrifuge database successfully built at {db_path}")


def _sanitize_seqid(seq_id: str) -> str:
    return "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in seq_id)[:80]


def kaiju_protein_id(taxid: str, accession: str = "") -> str:
    """Kaiju reads taxid as the integer after the last ``_`` in the FASTA id.

    Use ``{accession}_{taxid}``. A bare numeric id also works.
    """
    taxid = str(taxid).split(".")[0]
    acc = _sanitize_seqid(accession).strip("_") if accession else ""
    if not acc:
        return taxid
    return f"{acc}_{taxid}"


def _parse_gff_attributes(attr_field: str) -> dict:
    attrs = {}
    if "ID=" in attr_field or "Parent=" in attr_field or "protein_id=" in attr_field:
        for part in attr_field.strip().split(";"):
            if "=" in part:
                key, value = part.split("=", 1)
                attrs[key.strip()] = value.strip().strip('"')
        return attrs
    for part in attr_field.strip().split(";"):
        part = part.strip()
        if " " in part:
            key, value = part.split(" ", 1)
            attrs[key.strip()] = value.strip().strip('"')
    return attrs


def translate_cds_from_gff(fasta_path: str, gff_path: str, taxid: str, min_aa: int = 10) -> List[SeqRecord]:
    """Translate CDS features from GFF/GTF using a nucleotide FASTA."""
    sequences = SeqIO.to_dict(iter_seqio_fasta(fasta_path))
    cds_groups: Dict[str, list] = {}
    with open_text(gff_path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            seqid, _source, ftype, start, end, _score, strand, phase, attrs = fields[:9]
            if ftype.upper() not in {"CDS", "CDS_FEATURE"}:
                continue
            parsed = _parse_gff_attributes(attrs)
            protein_id = (
                parsed.get("protein_id")
                or parsed.get("ID")
                or parsed.get("Parent")
                or parsed.get("transcript_id")
                or f"{seqid}:{start}-{end}"
            )
            try:
                phase_i = int(phase) if phase not in {".", ""} else 0
            except ValueError:
                phase_i = 0
            cds_groups.setdefault(protein_id, []).append(
                (seqid, int(start), int(end), strand, phase_i)
            )

    records: List[SeqRecord] = []
    for protein_id, parts in cds_groups.items():
        strand = parts[0][3]
        ordered = sorted(parts, key=lambda x: x[1], reverse=(strand == "-"))
        cds = Seq("")
        for seqid, start, end, feat_strand, phase in ordered:
            record = sequences.get(seqid)
            if record is None:
                record = next((s for s in sequences.values() if seqid in s.id), None)
            if record is None:
                cds = Seq("")
                break
            chunk = record.seq[start - 1:end]
            if feat_strand == "-":
                chunk = chunk.reverse_complement()
            if phase:
                chunk = chunk[phase:]
            cds += chunk
        trim = (len(cds) // 3) * 3
        if trim < 3:
            continue
        protein = str(cds[:trim].translate(to_stop=True)).replace("*", "").replace("X", "")
        if len(protein) < min_aa:
            continue
        records.append(SeqRecord(Seq(protein), id=kaiju_protein_id(taxid, protein_id), description=""))
    return records


def nucleotide_to_frame_records(input_file: str, taxid: str) -> List[SeqRecord]:
    """One 6-frame translation per contig with stops stored as X.

    Kaiju-mkbwt treats ``*`` as the sequence terminator, so leftover stop
    symbols in the protein FASTA break the index. Headers must end with
    ``_{taxid}`` (Kaiju reads the integer after the last underscore).
    """
    records: List[SeqRecord] = []
    taxid = str(taxid).split(".")[0]
    for seq_index, nucl in enumerate(iter_seqio_fasta(input_file)):
        seq = str(nucl.seq).upper().replace("U", "T")
        strand_i = 0
        for strand_seq in (Seq(seq), Seq(seq).reverse_complement()):
            for frame in range(3):
                frame_seq = strand_seq[frame:]
                frame_seq = frame_seq[: (len(frame_seq) // 3) * 3]
                if len(frame_seq) < 3:
                    continue
                aa = str(frame_seq.translate()).replace("*", "X")
                if not aa:
                    continue
                records.append(
                    SeqRecord(
                        Seq(aa),
                        id=kaiju_protein_id(taxid, f"c{seq_index}s{strand_i}f{frame}"),
                        description="",
                    )
                )
            strand_i += 1
    return records


def nucleotide_to_orf_records(input_file: str, taxid: str, min_aa: int = 15) -> List[SeqRecord]:
    """Translate nucleotide FASTA in 6 frames (stop-free ORFs, taxid-only headers)."""
    records: List[SeqRecord] = []
    taxid = str(taxid).split(".")[0]
    for nucl in iter_seqio_fasta(input_file):
        seq = str(nucl.seq).upper().replace("U", "T")
        for strand_seq in (Seq(seq), Seq(seq).reverse_complement()):
            for frame in range(3):
                frame_seq = strand_seq[frame:]
                frame_seq = frame_seq[: (len(frame_seq) // 3) * 3]
                if len(frame_seq) < 3:
                    continue
                aa = str(frame_seq.translate())
                for orf in aa.split("*"):
                    orf = orf.replace("X", "")
                    if len(orf) < min_aa:
                        continue
                    records.append(
                        SeqRecord(Seq(orf), id=kaiju_protein_id(taxid, f"orf{len(records)}"), description="")
                    )
    return records


def find_local_proteome_or_gff(input_file: str, taxid: str) -> Tuple[Optional[str], Optional[str]]:
    """Return ('protein'|'gff', path) for a sibling annotation file if present."""
    for cand in sibling_candidates(input_file, taxid, (".faa", ".faa.gz")):
        if cand.exists() and (cand.suffix.lower() == ".faa" or cand.name.lower().endswith(".faa.gz")):
            return "protein", str(cand)
    for cand in sibling_candidates(
        input_file, taxid, (".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz")
    ):
        if cand.exists():
            return "gff", str(cand)
    return None, None


def process_fasta_kaiju(input_file: str, 
                       taxid: str,
                       db_path: str = "kaiju_db",
                       protein: bool = False,
                       email: Optional[str] = None,
                       fetch_missing: bool = True,
                       min_aa: int = 30) -> str:
    """Process FASTA into a Kaiju protein library.

    Kaiju indexes proteins and reads taxid as the integer after the last ``_``
    in the FASTA header (``{accession}_{taxid}``). Nucleotide genomes become:
      1. sibling `.faa` proteome, or
      2. sibling GFF/GTF CDS translation, or
      3. downloaded NCBI proteome / GFF, or
      4. 6-frame translations of the provided FASTA.
    """
    email = email or default_entrez_email()
    taxid = str(taxid).split(".")[0]
    temp_fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.faa', delete=False)
    temp_fasta.close()
    records: List[SeqRecord] = []

    suffix = Path(str(input_file).removesuffix(".gz")).suffix.lower()
    is_protein = protein or suffix == ".faa"
    if is_protein:
        for rec in iter_seqio_fasta(input_file):
            aa = str(rec.seq).replace("*", "").replace("X", "")
            if len(aa) < min_aa:
                continue
            records.append(SeqRecord(Seq(aa), id=kaiju_protein_id(taxid, rec.id), description=""))
    else:
        kind, annot_path = find_local_proteome_or_gff(input_file, taxid)
        if kind is None and fetch_missing:
            from samovar.genome_index import proteome_storage_dir, resolve_indexed_path

            cache_dir = str(proteome_storage_dir())
            Path(cache_dir).mkdir(parents=True, exist_ok=True)
            indexed = resolve_indexed_path(taxid)
            if indexed is not None:
                sibling_kind, sibling_path = find_local_proteome_or_gff(str(indexed), taxid)
                if sibling_kind is not None:
                    kind, annot_path = sibling_kind, sibling_path
            if kind is None:
                proteome = fetch_proteome(taxid, cache_dir, email, silent=False)
                if proteome:
                    kind, annot_path = "protein", proteome
                else:
                    gff = fetch_gff(taxid, cache_dir, email, silent=False)
                    if gff:
                        kind, annot_path = "gff", gff

        if kind == "protein":
            for rec in iter_seqio_fasta(annot_path):
                aa = str(rec.seq).replace("*", "").replace("X", "")
                if len(aa) < min_aa:
                    continue
                records.append(SeqRecord(Seq(aa), id=kaiju_protein_id(taxid, rec.id), description=""))
        elif kind == "gff":
            records.extend(translate_cds_from_gff(input_file, annot_path, taxid, min_aa=min_aa))

        # Index 6-frame translations of the provided nucleotide file so reads
        # simulated from that FASTA can still be classified.
        records.extend(nucleotide_to_frame_records(input_file, taxid))
        if kind is None:
            logger.warning(
                "No proteome/GFF CDS found for %s (taxid %s); using 6-frame translations",
                input_file,
                taxid,
            )

    if not records:
        raise ValueError(f"No protein sequences produced for {input_file} (taxid {taxid})")

    SeqIO.write(records, temp_fasta.name, "fasta")
    return temp_fasta.name

def add_database_kaiju(
    input_file: str,
    taxid: str,
    db_path: str = "kaiju_db",
    protein: bool = False,
    email: Optional[str] = None,
    fetch_missing: bool = True,
) -> None:
    """Add sequences to the Kaiju database library.
    
    This function processes a FASTA file and adds its sequences to the Kaiju database library
    with proper taxonomy IDs. The sequences will be translated to protein if they are not already.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to associate with the sequences.
        db_path (str, optional): Path to the Kaiju database. Defaults to "kaiju_db".
        protein (bool, optional): Whether the input is already protein sequences. Defaults to False.
    
    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    # Create database directory if it doesn't exist
    os.makedirs(db_path, exist_ok=True)
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Process the input file
    processed_file = process_fasta_kaiju(
        input_file, taxid, db_path, protein, email=email, fetch_missing=fetch_missing
    )
    
    # Append to the library file (gunzip a previous compressed library if needed)
    library_file = os.path.join(db_path, "library.faa")
    library_gz = library_file + ".gz"
    if not os.path.exists(library_file) and os.path.exists(library_gz):
        gunzip_file(library_gz, library_file, remove_source=True)
    with open(processed_file, 'r') as infile, open(library_file, 'a') as outfile:
        outfile.write(infile.read())
    
    # Clean up temporary file
    os.remove(processed_file)
    
    logger.info(f"Sequences successfully added to Kaiju database at {db_path}")

def build_database_kaiju(
    db_path: str = "kaiju_db",
    threads: int = 1,
    protein: bool = False
) -> None:
    """Build a custom Kaiju database from the processed sequences.
    
    This function builds a Kaiju database using the sequences that have been processed
    and added to the database directory. It builds the BWT first, then the FM-index
    for efficient sequence classification.
    
    Args:
        db_path (str, optional): Path to the database directory. Defaults to "kaiju_db".
        threads (int, optional): Number of threads to use for building. Defaults to 1.
        protein (bool, optional): Whether the input is already protein sequences. Defaults to False.
    
    Note:
        Building a large database can be memory-intensive and time-consuming.
    """
    # Create database directory if it doesn't exist
    os.makedirs(db_path, exist_ok=True)
    
    library_file = os.path.join(db_path, "library.faa")
    library_gz = library_file + ".gz"
    if not os.path.exists(library_file) and os.path.exists(library_gz):
        gunzip_file(library_gz, library_file, remove_source=False)
    base_name = os.path.join(db_path, "kaiju_db")
    
    # First build the BWT
    mkbwt_cmd = [
        "kaiju-mkbwt",
        "-n", str(threads),
        "-a", "ACDEFGHIKLMNPQRSTVWYX",
        "-o", base_name,
        library_file,
    ]
    run_command(mkbwt_cmd)

    bwt_file = f"{base_name}.bwt"
    if not os.path.exists(bwt_file) and os.path.exists(f"{base_name}.fmi"):
        logger.info("kaiju-mkbwt already produced an FM-index; skipping kaiju-mkfmi")
        if os.path.exists(library_file):
            gzip_file(library_file)
        return
    if not os.path.exists(bwt_file):
        raise FileNotFoundError(
            f"kaiju-mkbwt did not write {bwt_file}. Check kaiju-mkbwt output and disk space."
        )

    mkfmi_cmd = [
        "kaiju-mkfmi",
        base_name
    ]
    run_command(mkfmi_cmd)
    
    # Clean up temporary files
    sa_file = f"{base_name}.sa"
    if os.path.exists(sa_file):
        os.remove(sa_file)
    
    bwt_file = f"{base_name}.bwt"
    if os.path.exists(bwt_file):
        os.remove(bwt_file)

    # Classification uses the FM-index, not the protein library.
    if os.path.exists(library_file):
        gzip_file(library_file)
    
    logger.info(f"Kaiju database successfully built at {db_path}")

class KrakenUniqLibrarySkip(Exception):
    """One taxon cannot be indexed. The rest of the library can still be built."""

    def __init__(self, taxid: str, path: str, reason: str):
        self.taxid = str(taxid)
        self.path = str(path)
        self.reason = str(reason)
        super().__init__(f"{self.taxid}\t{self.reason}\t{self.path}")


def _krakenuniq_taxid_name(taxid: str) -> str:
    text = str(taxid or "").strip().replace("/", "_").replace("\\", "_")
    return text or "unknown"


def _is_protein_fasta(path: str) -> bool:
    name = str(path or "").lower()
    return name.endswith(".faa") or name.endswith(".faa.gz")


def krakenuniq_added_dir(db_path: str) -> Path:
    """``library/added/`` — one uncompressed FASTA per taxid, as krakenuniq-build scans it."""
    added = Path(db_path) / "library" / "added"
    added.mkdir(parents=True, exist_ok=True)
    return added


def process_fasta_krakenunique(input_file: str, taxid: str, genome_name: str, db_path: str) -> tuple[str, str]:
    """Write one uncompressed nucleotide FASTA and its seqid map under ``library/added/``.

    ``krakenuniq-build`` finds ``*.fa`` under ``library/`` and concatenates
    ``library/**/*.map`` into ``seqid2taxid.map``. Protein FASTA and missing
    files are not invented; they raise ``KrakenUniqLibrarySkip``.

    Returns:
        Paths to ``library/added/{taxid}.fa`` and ``library/added/{taxid}.map``.
    """
    if _is_protein_fasta(input_file):
        raise KrakenUniqLibrarySkip(taxid, input_file, "proteome")
    if not os.path.exists(input_file):
        raise KrakenUniqLibrarySkip(taxid, input_file, "missing")

    added = krakenuniq_added_dir(db_path)
    stem = _krakenuniq_taxid_name(taxid)
    processed_fasta = str(added / f"{stem}.fa")
    preprocess_fasta(input_file, processed_fasta, mutation_rate=0, include_percent=100)

    map_file = str(added / f"{stem}.map")
    wrote = 0
    with open(map_file, "w", encoding="utf-8") as map_out:
        for record in iter_seqio_fasta(processed_fasta):
            seq_id = record.id.split()[0]
            if not seq_id or len(record.seq) == 0:
                continue
            map_out.write(f"{seq_id}\t{taxid}\t{genome_name}\n")
            wrote += 1
    if wrote == 0:
        for path in (processed_fasta, map_file):
            try:
                os.remove(path)
            except OSError:
                pass
        raise KrakenUniqLibrarySkip(taxid, input_file, "empty")
    return processed_fasta, map_file

def add_database_krakenunique(
    input_file: str,
    taxid: str,
    genome_name: str,
    db_path: str = "krakenuniq_db"
) -> tuple[str, str]:
    """Add one nucleotide genome to ``library/added/{taxid}.fa``.

    Raises:
        KrakenUniqLibrarySkip: Genome file is missing, empty, or protein-only.
    """
    os.makedirs(db_path, exist_ok=True)
    processed, map_file = process_fasta_krakenunique(input_file, taxid, genome_name, db_path)
    logger.info("Added KrakenUniq genome %s (%s) to %s", taxid, genome_name, processed)
    return processed, map_file


def write_krakenuniq_skips(db_path: str, skips: List[KrakenUniqLibrarySkip]) -> Path:
    """Taxa left out of the library. Not placed under ``library/``, so the build does not index them."""
    dest = Path(db_path) / "missing-taxa.tsv"
    dest.parent.mkdir(parents=True, exist_ok=True)
    lines = ["taxid\treason\tpath"]
    for item in skips:
        lines.append(f"{item.taxid}\t{item.reason}\t{item.path}")
    dest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return dest

def build_database_krakenunique(
    db_path: str = "krakenuniq_db",
    threads: int = 1,
    minimizer_len: Optional[int] = None,
) -> None:
    """Build a custom KrakenUniq database from the library of sequences.
    
    This function builds a KrakenUniq database using the sequences that have been added to the library.
    
    Args:
        db_path (str, optional): Path to the database directory. Defaults to "krakenuniq_db".
        threads (int, optional): Number of threads to use for building. Defaults to 1.
        minimizer_len: Override the library-sized minimizer index (12–15 by default).
    
    Note:
        Building a large database can be memory-intensive and time-consuming.
    """
    from samovar.paths import resolve_executable

    def _on_disk(name: str) -> str:
        head = command_argv(resolve_executable(name) or "")
        path = Path(head[0]) if head else Path()
        if path.is_file() and os.access(path, os.X_OK):
            return str(path)
        return ""

    build_bin = _on_disk("krakenuniq-build")
    if not build_bin:
        raise RuntimeError(
            "krakenuniq-build not found. Import it next to krakenuniq (./install.sh full ask)."
        )
    jellyfish_path = _on_disk("jellyfish")
    if not jellyfish_path:
        raise RuntimeError("jellyfish not found. Import it next to krakenuniq (./install.sh full ask).")

    library = Path(db_path) / "library"
    library_files = list(library.rglob("*.fa")) + list(library.rglob("*.fna"))
    if not library_files:
        raise FileNotFoundError(
            f"No nucleotide FASTA under {library}. "
            "KrakenUniq indexes library/*.fa; see missing-taxa.tsv for genomes that were skipped."
        )
    if minimizer_len is None:
        library_bytes = sum(path.stat().st_size for path in library_files)
        minimizer_len = 12
        # db_sort allocates 8 * 4**length bytes even for a tiny library.
        # Grow beyond the 128 MiB baseline only when the library is as large.
        for length in range(13, 16):
            if 8 * 4**length > library_bytes:
                break
            minimizer_len = length
    if not 1 <= minimizer_len < 31:
        raise ValueError("KrakenUniq minimizer length must be between 1 and 30")
    taxonomy = Path(db_path) / "taxonomy"
    if not (taxonomy / "nodes.dmp").is_file() or not (taxonomy / "names.dmp").is_file():
        raise FileNotFoundError(
            f"KrakenUniq needs {taxonomy}/nodes.dmp and names.dmp "
            "(copied from the configured taxdump)."
        )
    # An empty or stale list hides the library. build_db.sh rebuilds it when absent.
    stale_list = Path(db_path) / "library-files.txt"
    if stale_list.exists():
        stale_list.unlink()

    # Set up environment variables
    env = os.environ.copy()
    env['JELLYFISH_BIN'] = jellyfish_path

    build_cmd = [
        build_bin,
        # Existing indexes/maps may contain synthetic IDs or a different index
        # length. Rebuild from library/ rather than silently reusing those files.
        "--rebuild",
        "--db", db_path,
        "--threads", str(threads),
        "--minimizer-len", str(minimizer_len),
    ]
    
    # Run command with modified environment
    logger.info(f"Running command: {' '.join(build_cmd)}")
    result = subprocess.run(build_cmd, check=True, text=True, env=env)
    
    logger.info(f"KrakenUniq database successfully built at {db_path}")

def build_database_from_config(
    config_path: str,
    db_type: str = "kaiju",
    db_path: str = "tests_outs/db",
    example_omit: Optional[bool] = None,
    threads: int = 1,
    minimizer_len: Optional[int] = None,
):
    """
    Build database from config file.
    
    Args:
        config_path: Path to config YAML file
        db_type: Type of database to build ("kaiju", "kraken2", "krakenunique", or "centrifuge")
        db_path: Path to store the database
        example_omit: If True, apply toy-only taxon gaps (Escherichia off Kraken2,
            Phage Phi X off Kaiju). If None, honour SAMOVAR_EXAMPLE_OMIT or
            auto-detect ``test_genomes`` input directories.
        threads: Workers for ``kaiju-mkbwt -n``, ``kraken2-build --threads``,
            and ``centrifuge-build -p``. Same count as ``samovar prepare --cores``.
        minimizer_len: Optional KrakenUniq minimizer length; default scales with the library.
    """
    threads = max(1, int(threads))
    from samovar.example_db import (
        filter_example_omit,
        should_apply_example_omit,
        warn_example_omit,
    )

    # Load configuration
    with open(config_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # Process input directories
    input_dirs = config['input_dir']
    if isinstance(input_dirs, str):
        input_dirs = [input_dirs]
    file_taxid_map = process_fasta_directories(input_dirs)
    from samovar.genome_index import is_assembly_accession, numeric_taxid_for
    from samovar.genome_fetcher import assembly_taxid, default_entrez_email

    resolved_map = {}
    email = default_entrez_email()
    for path, taxid in file_taxid_map.items():
        ident = str(taxid)
        resolved = numeric_taxid_for(ident)
        if is_assembly_accession(ident) and not str(resolved).split(".")[0].isdigit():
            looked = assembly_taxid(ident, email)
            if looked:
                resolved = looked
        resolved_map[path] = resolved
    file_taxid_map = resolved_map
    preferred = {}
    for path, taxid in file_taxid_map.items():
        prev = preferred.get(taxid)
        if prev is None:
            preferred[taxid] = path
            continue
        # Prefer nucleotide genomes so proteome + ORF indexing can both run.
        nucleotide = (
            ".fna", ".fa", ".fasta", ".fsa", ".ffn",
            ".fna.gz", ".fa.gz", ".fasta.gz", ".fsa.gz", ".ffn.gz",
        )
        if path.endswith(nucleotide) and not prev.endswith(nucleotide):
            preferred[taxid] = path
    file_taxid_map = {path: taxid for taxid, path in preferred.items()}

    apply_omit = should_apply_example_omit(input_dirs, explicit=example_omit)
    omitted: Dict[str, str] = {}
    if apply_omit and db_type in {"kaiju", "kraken2"}:
        file_taxid_map, omitted = filter_example_omit(file_taxid_map, db_type)
        warn_example_omit(db_type, omitted.values())
    
    # Get lists of files and taxids
    input_files = list(file_taxid_map.keys())
    taxids = list(file_taxid_map.values())
    
    print(f"Processing {len(input_files)} files with taxids: {taxids}")
    if omitted:
        print(f"EXAMPLE omit ({db_type}): skipped taxids {sorted(set(omitted.values()))}")
    if not input_files:
        raise ValueError(
            f"No genomes left to index for {db_type} after example omit. "
            "This flag is only for toy tests; use --no-example-omit for a full index."
        )

    
    # Build database based on type
    if db_type == "kaiju":
        library_file = os.path.join(db_path, "library.faa")
        if os.path.exists(library_file):
            os.remove(library_file)
        library_gz = library_file + ".gz"
        if os.path.exists(library_gz):
            os.remove(library_gz)
        for input_file, taxid in zip(input_files, taxids):
            add_database_kaiju(
                input_file,
                taxid,
                db_path=db_path,
                protein=str(input_file).endswith(".faa") or str(input_file).endswith(".faa.gz"),
                fetch_missing=True,
            )
        get_taxonomy_db(db_path=db_path)
        build_database_kaiju(db_path=db_path, threads=threads, protein=False)
    elif db_type == "kraken2":
        work_dir, final_dir = kraken2_db_workdir(db_path)
        try:
            for input_file, taxid in zip(input_files, taxids):
                add_database_kraken2(input_file, taxid, db_path=work_dir)
            get_taxonomy_db(db_path=work_dir + "/taxonomy")
            build_database_kraken2(
                db_path=work_dir,
                threads=threads,
                kmer_len=35,
                minimizer_len=31,
                minimizer_spaces=7,
                skip_maps=True,
            )
            if final_dir is not None:
                install_kraken2_db(work_dir, final_dir)
                work_dir = ""
        finally:
            if work_dir and final_dir is not None and os.path.isdir(work_dir):
                shutil.rmtree(work_dir, ignore_errors=True)
    elif db_type == "krakenunique":
        os.makedirs(db_path, exist_ok=True)
        skips: List[KrakenUniqLibrarySkip] = []
        added = 0
        for input_file, taxid in zip(input_files, taxids):
            genome_name = os.path.basename(input_file)
            try:
                add_database_krakenunique(input_file, taxid, genome_name, db_path=db_path)
            except KrakenUniqLibrarySkip as exc:
                skips.append(exc)
                logger.warning(
                    "KrakenUniq skipping taxid %s (%s): %s",
                    exc.taxid,
                    exc.reason,
                    exc.path,
                )
                continue
            added += 1
        write_krakenuniq_skips(db_path, skips)
        if added == 0:
            raise ValueError(
                "No nucleotide genomes left to index for krakenunique. "
                f"Skipped {len(skips)} taxa; see {db_path}/missing-taxa.tsv."
            )
        # Dumps belong in taxonomy/, not next to the FASTA. Reuse the configured taxdump.
        get_taxonomy_db(db_path=os.path.join(db_path, "taxonomy"))
        build_database_krakenunique(db_path=db_path, threads=threads, minimizer_len=minimizer_len)
    elif db_type == "centrifuge":
        os.makedirs(db_path, exist_ok=True)
        nucleotide = [
            (path, taxid)
            for path, taxid in zip(input_files, taxids)
            if _is_nucleotide_genome(path)
        ]
        if not nucleotide:
            raise ValueError("No nucleotide genomes left to index for centrifuge.")
        written = write_centrifuge_inputs(
            nucleotide,
            os.path.join(db_path, "library.fna"),
            os.path.join(db_path, "seqid2taxid.map"),
        )
        if written == 0:
            raise ValueError("No nucleotide sequences left to index for centrifuge.")
        get_taxonomy_db(db_path=os.path.join(db_path, "taxonomy"))
        build_database_centrifuge(db_path=db_path, threads=threads)
    elif db_type in {"kmer_encoder", "metauto", "kmerenc"}:
        from samovar.kmer_encoder import train_encoder

        os.makedirs(db_path, exist_ok=True)
        dest = os.path.join(db_path, "kmer_encoder.pt")
        train_encoder(input_dirs, dest)
        logger.info("k-mer encoder model written to %s", dest)
    else:
        raise ValueError(f"Unknown database type: {db_type}")
