#!/usr/bin/env python3

import os
import subprocess
import logging
import tempfile
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
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
    logger.info(f"Running command: {' '.join(cmd)}")
    return subprocess.run(cmd, check=check, text=True)

def get_taxonomy_db(
    db_path: str="kraken_db",
    taxonomy_path: str=None
) -> None:
    """Download or copy the Kraken2 taxonomy database.
    
    This function either downloads the NCBI taxonomy database or copies it from a specified path.
    The taxonomy database is required for Kraken2 to function properly.
    
    Args:
        db_path (str, optional): Path where the database will be stored. Defaults to "kraken_db".
        taxonomy_path (str, optional): Path to existing taxonomy database. If None, downloads from NCBI. Defaults to None.
    
    Note:
        If taxonomy_path is None, the function will download the latest taxonomy database from NCBI.
        The downloaded file will be automatically extracted to the specified db_path.
        Override the URL with SAMOVAR_TAXDUMP_URL.
    """
    os.makedirs(db_path, exist_ok=True)
    existing_nodes = Path(db_path) / "nodes.dmp"
    existing_nested = Path(db_path) / "taxonomy" / "nodes.dmp"
    if existing_nodes.exists() or existing_nested.exists():
        logger.info(f"Taxonomy already present under {db_path}; skipping download")
        return

    if taxonomy_path is None:
        taxdump_url = os.environ.get(
            "SAMOVAR_TAXDUMP_URL",
            "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz",
        )
        download_cmd =[
            "wget",
            taxdump_url,
            "-P",
            db_path
        ]
        run_command(download_cmd)

        extract_cmd = [
            "tar",
            "-xzf",
            db_path + "/taxdump.tar.gz",
            "-C",
            db_path
        ]
        run_command(extract_cmd)

        logger.info(f"Taxonomy database downloaded and extracted to {db_path}")
    
    else:
        copy_cmd = [
            "cp",
            "-r",
            taxonomy_path,
            db_path
        ]
        run_command(copy_cmd)

        logger.info(f"Taxonomy database copied to {db_path}/taxonomy")

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
            cache_dir = str(Path(db_path) / "proteomes")
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

def process_fasta_krakenunique(input_file: str, taxid: str, genome_name: str, db_path: str) -> tuple[str, str]:
    """Process FASTA file for KrakenUniq database building.
    
    This function preprocesses a FASTA file and creates a genomes.map file entry.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to be added to sequence headers.
        genome_name (str): Name of the genome/assembly.
        db_path (str): Path to the database directory.
    
    Returns:
        tuple[str, str]: Paths to the processed FASTA file and genomes.map file.
    """
    # Preprocess the FASTA file
    processed_fasta = os.path.join(db_path, os.path.basename(input_file))
    preprocess_fasta(input_file, processed_fasta, mutation_rate=0, include_percent=100)
    
    # Create genomes.map entry
    map_file = os.path.join(db_path, "genomes.map")
    with open(map_file, 'a') as map_out:
        for record in iter_seqio_fasta(processed_fasta):
            # Get sequence ID (header without '>' up to first space)
            seq_id = record.id.split()[0]
            # Write to genomes.map
            map_out.write(f"{seq_id}\t{taxid}\t{genome_name}\n")
    
    return processed_fasta, map_file

def add_database_krakenunique(
    input_file: str,
    taxid: str,
    genome_name: str,
    db_path: str = "krakenuniq_db"
) -> None:
    """Add sequences to the KrakenUniq database library.
    
    This function processes a FASTA file and adds its sequences to the KrakenUniq database library
    with proper taxonomy IDs and creates a genomes.map file.
    
    Args:
        input_file (str): Path to the input FASTA file.
        taxid (str): Taxonomy ID to associate with the sequences.
        genome_name (str): Name of the genome/assembly.
        db_path (str, optional): Path to the KrakenUniq database. Defaults to "krakenuniq_db".
    
    Raises:
        FileNotFoundError: If the input file does not exist.
    """
    # Create database directory if it doesn't exist
    os.makedirs(db_path, exist_ok=True)
    
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file not found: {input_file}")
    
    # Process the input file
    process_fasta_krakenunique(input_file, taxid, genome_name, db_path)
    
    logger.info(f"Sequences successfully added to KrakenUniq database at {db_path}")

def build_database_krakenunique(
    db_path: str = "krakenuniq_db",
    threads: int = 1
) -> None:
    """Build a custom KrakenUniq database from the library of sequences.
    
    This function builds a KrakenUniq database using the sequences that have been added to the library.
    
    Args:
        db_path (str, optional): Path to the database directory. Defaults to "krakenuniq_db".
        threads (int, optional): Number of threads to use for building. Defaults to 1.
    
    Note:
        Building a large database can be memory-intensive and time-consuming.
    """
    # Check if krakenuniq is installed
    try:
        subprocess.run(['which', 'krakenuniq-build'], check=True, capture_output=True)
    except subprocess.CalledProcessError:
        raise RuntimeError("krakenuniq-build not found. Please install KrakenUniq first.")

    # Check if jellyfish is installed
    try:
        jellyfish_path = subprocess.run(['which', 'jellyfish'], check=True, capture_output=True, text=True).stdout.strip()
    except subprocess.CalledProcessError:
        raise RuntimeError("jellyfish not found. Please install Jellyfish first.")

    # Set up environment variables
    env = os.environ.copy()
    env['JELLYFISH_BIN'] = jellyfish_path

    build_cmd = [
        "krakenuniq-build",
        "--db", db_path,
        "--threads", str(threads),
        "--taxids-for-genomes",
        "--taxids-for-sequences"
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
):
    """
    Build database from config file.
    
    Args:
        config_path: Path to config YAML file
        db_type: Type of database to build ("kaiju", "kraken2", or "krakenunique")
        db_path: Path to store the database
        example_omit: If True, apply toy-only taxon gaps (Escherichia off Kraken2,
            Phage Phi X off Kaiju). If None, honour SAMOVAR_EXAMPLE_OMIT or
            auto-detect ``test_genomes`` input directories.
    """
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
        build_database_kaiju(db_path=db_path, threads=1, protein=False)
    elif db_type == "kraken2":
        for input_file, taxid in zip(input_files, taxids):
            add_database_kraken2(input_file, taxid, db_path=db_path)
        get_taxonomy_db(db_path=db_path+"/taxonomy")
        build_database_kraken2(db_path=db_path, threads=1, kmer_len=35, 
                             minimizer_len=31, minimizer_spaces=7, skip_maps=True)
    elif db_type == "krakenunique":
        # Create output directory
        os.makedirs(db_path, exist_ok=True)
        # Get taxonomy database
        get_taxonomy_db(db_path=db_path)
        # Add each genome to the database
        for input_file, taxid in zip(input_files, taxids):
            # Use filename as genome name
            genome_name = os.path.basename(input_file)
            add_database_krakenunique(input_file, taxid, genome_name, db_path=db_path)
        # Build the database
        build_database_krakenunique(db_path=db_path, threads=1)
    else:
        raise ValueError(f"Unknown database type: {db_type}")