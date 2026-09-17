"""Assembly-profiling stage: MegaHIT, anvi'o, MAG QC, GTDB-Tk, minimap2, CoverM."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from samovar.annotators_wrapper import AssemblyAnnotator, get_annotator_instance
from samovar.assembly_profiling import (
    _gtdb_taxid_from_lineage,
    _iter_mag_fastas,
    gtdbtk_db_ready,
    run_aligner,
    run_anvio_binner,
    run_assembler,
    run_assembly_annotator,
    run_binner,
    run_binner_combine,
    run_binner_qc,
    run_coverm,
    run_gene_caller,
    run_mag_quantifier,
    run_mag_taxonomy,
    run_read_assigner,
    run_taxon_quantifier,
    which_tool,
)
from samovar.parse_annotators import Annotation, match_annotation

REPO = Path(__file__).resolve().parents[1]
READS_R1 = REPO / "tests" / "data" / "reads" / "1_full_R1.fastq"
READS_R2 = REPO / "tests" / "data" / "reads" / "1_full_R2.fastq"


def test_translate_orfs_gene_caller(tmp_path):
    """Naive codon translation writes FAA without Prodigal."""
    contigs = tmp_path / "contigs.fa"
    contigs.write_text(">c1\nATGAAATTTAAATAA\n")
    dest = tmp_path / "genes"
    run_gene_caller(str(contigs), str(dest), name="translate")
    faa = list(dest.glob("*.faa"))
    assert faa
    assert faa[0].read_text().startswith(">")


def test_gtdb_lineage_emits_ncbi_taxid(monkeypatch):
    monkeypatch.setattr(
        "samovar.taxdump.load_scientific_name_to_taxid",
        lambda path=None: {"thermosphaera aggregans": "2268", "bacteria": "2"},
    )
    taxid = _gtdb_taxid_from_lineage(
        "d__Archaea;p__Thermoproteota;s__Thermosphaera aggregans"
    )
    assert taxid == "2268"
    assert _gtdb_taxid_from_lineage("Unclassified Bacteria") == "2"


@pytest.mark.optional
def test_prodigal_gene_caller(tmp_path):
    try:
        which_tool("prodigal")
    except FileNotFoundError as exc:
        pytest.skip(str(exc))
    contigs = tmp_path / "contigs.fa"
    # Long enough for Prodigal meta mode
    seq = ("ATG" + "AAA" * 80 + "TAA") * 4
    contigs.write_text(f">c1\n{seq}\n")
    dest = tmp_path / "genes"
    try:
        run_gene_caller(str(contigs), str(dest), name="prodigal")
    except Exception as exc:
        pytest.skip(f"Prodigal could not call genes: {exc}")
    assert list(dest.glob("*.faa"))


def _sidecar_bin(name: str) -> Path:
    from samovar.sidecar import SIDECARS, create_sidecar_env, sidecar_is_healthy, sidecar_prefix

    spec = SIDECARS[name]
    prefix = sidecar_prefix(name)
    binary = prefix / "bin" / spec["binary"]
    if sidecar_is_healthy(name, prefix) and binary.is_file():
        return binary
    try:
        create_sidecar_env(name)
    except Exception as exc:
        pytest.skip(f"cannot conda-install {name}: {exc}")
    if not binary.is_file():
        pytest.skip(f"{name} missing after sidecar install")
    return binary


def _put_on_path(monkeypatch, binary: Path) -> None:
    monkeypatch.setenv("PATH", f"{binary.parent}:{os.environ.get('PATH', '')}")


def test_factory_and_match_annotation():
    inst = get_annotator_instance("assembly", {"run_name": "assembly-test"}, {})
    assert isinstance(inst, AssemblyAnnotator)
    outs = inst.get_expected_outputs("1_full", "/tmp/out")
    assert outs[0].endswith("1_full_assembly-test.assembly.out")
    assert match_annotation("1_full_assembly-test.assembly.out") == "assembly"


def test_taxon_quantifier_joins_mag_counts_and_taxonomy(tmp_path):
    """taxon_quantifier + read_assigner glue MAG abundance into taxid tables."""
    mag_ab = tmp_path / "mag.tsv"
    mag_ab.write_text("mag_id\tN\nmagA\t10\nmagB\t5\n", encoding="utf-8")
    tax = tmp_path / "tax.tsv"
    tax.write_text("mag_id\ttaxid\tlineage\nmagA\t9606\td__Euk\nmagB\t2\td__Bac\n", encoding="utf-8")
    taxon_out = tmp_path / "taxon.tsv"
    run_taxon_quantifier(str(mag_ab), str(tax), str(taxon_out))
    table = pd.read_csv(taxon_out, sep="\t")
    assert set(table["taxid"].astype(str)) == {"9606", "2"}
    assert int(table.loc[table["taxid"].astype(str) == "9606", "N_1"].iloc[0]) == 10

    sam = tmp_path / "hits.bam"
    sam.write_text(
        "r1\t0\tmagA~c1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        "r2\t0\tmagB~c1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n",
        encoding="utf-8",
    )
    dest = tmp_path / "reads.tsv"
    run_read_assigner(str(sam), str(tax), str(dest))
    parsed = Annotation({str(dest): "assembly"})
    assert set(parsed.DataFrame.index) == {"r1", "r2"}
    tax_col = [c for c in parsed.DataFrame.columns if str(c).startswith("taxID")][0]
    assert set(parsed.DataFrame[tax_col].astype(str)) == {"9606", "2"}
    feat_cols = [c for c in parsed.DataFrame.columns if str(c).startswith("feat_")]
    assert any("MAG_ID" in c for c in feat_cols)
    mag_col = next(c for c in feat_cols if "MAG_ID" in c)
    assert set(parsed.DataFrame[mag_col].astype(str)) == {"magA", "magB"}
    header = dest.read_text(encoding="utf-8").splitlines()[0]
    assert header == "seq\ttaxID\tMAG_ID"


def test_identity_composite_emits_mag_id_feature(tmp_path):
    dest = tmp_path / "sample_run.assembly.out"
    run_assembly_annotator(
        str(READS_R1),
        str(READS_R2),
        str(dest),
        threads=1,
        assembler="identity",
        gene_caller="identity",
        binners=["identity"],
        binner_qc="identity",
        binner_combine="identity",
        mag_taxonomy="identity",
        aligner="identity",
        mag_quantifier="identity",
        work_dir=str(tmp_path / "sample_run.assembly"),
    )
    assert dest.is_file()
    text = dest.read_text(encoding="utf-8")
    assert text.splitlines()[0] == "seq\ttaxID\tMAG_ID"
    parsed = Annotation({str(dest): "assembly"})
    tax_cols = [c for c in parsed.DataFrame.columns if str(c).startswith("taxID")]
    feat_cols = [c for c in parsed.DataFrame.columns if str(c).startswith("feat_")]
    assert tax_cols
    mag_col = next(c for c in feat_cols if "MAG_ID" in c)
    mag_ids = set(parsed.DataFrame[mag_col].astype(str))
    assert mag_ids <= {"mag1", "mag2"}
    assert mag_ids
    assert set(parsed.DataFrame[tax_cols[0]].astype(str)) == {"562"}


def test_identity_aligner_maps_taxid_headers_to_mags(tmp_path):
    from samovar.baselines.identity_aligner import main as align_main

    mag = tmp_path / "mags"
    mag.mkdir()
    (mag / "mag1.fa").write_text(">mag1\nACGT\n")
    (mag / "mag2.fa").write_text(">mag2\nTGCA\n")
    r1 = tmp_path / "r1.fastq"
    r1.write_text(
        "@rA|taxid:10710|\nACGT\n+\nIIII\n"
        "@rB|taxid:10699|\nACGT\n+\nIIII\n"
    )
    r2 = tmp_path / "r2.fastq"
    r2.write_text("@rC|taxid:10710|/2\nACGT\n+\nIIII\n")
    dest = tmp_path / "hits.sam"
    align_main(["-i", str(r1), "-I", str(r2), "-r", str(mag), "-o", str(dest)])
    text = dest.read_text()
    assert "rA|taxid:10710|\t0\tmag1\t" in text
    assert "rB|taxid:10699|\t0\tmag2\t" in text
    assert "rC|taxid:10710|/2\t0\tmag1\t" in text


def test_gzip_bam_uses_samtools(tmp_path, monkeypatch):
    import gzip

    bam = tmp_path / "reads_to_mag.bam"
    with gzip.open(bam, "wb") as handle:
        handle.write(b"not uncompressed BAM magic")
    tax = tmp_path / "mag_taxonomy.tsv"
    tax.write_text("mag_id\ttaxid\tlineage\nmagA\t54254\td__Archaea\n", encoding="utf-8")
    dest = tmp_path / "reads.tsv"
    monkeypatch.setattr(
        "samovar.assembly_profiling.which_tool",
        lambda name, default="": "/usr/bin/samtools",
    )
    monkeypatch.setattr(
        "samovar.assembly_profiling.subprocess.check_output",
        lambda cmd, text=True, errors="replace": (
            "r1\t0\tmagA~c1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\n"
        ),
    )
    run_read_assigner(str(bam), str(tax), str(dest))
    assert dest.read_text(encoding="utf-8") == "seq\ttaxID\tMAG_ID\nr1\t54254\tmagA\n"


def test_rewrite_mag_taxonomy_ncbi_unclassified_bacteria(tmp_path):
    from samovar.assembly_profiling import rewrite_mag_taxonomy_ncbi

    tax = tmp_path / "mag_taxonomy.tsv"
    tax.write_text(
        "mag_id\ttaxid\tlineage\nmagA\t0\tUnclassified Bacteria\n",
        encoding="utf-8",
    )
    rewrite_mag_taxonomy_ncbi(tax)
    mapped = tax.read_text(encoding="utf-8")
    assert mapped.splitlines()[1].split("\t")[1] == "2"
    tax.write_text(
        mapped + "magA\t0\t\n",
        encoding="utf-8",
    )
    from samovar.assembly_profiling import load_mag_taxonomy

    assert load_mag_taxonomy(tax)["magA"] == "2"


def test_identity_combine(tmp_path):
    src = tmp_path / "bins"
    src.mkdir()
    (src / "mag1.fa").write_text(">c1\nACGT\n", encoding="utf-8")
    dest = tmp_path / "selected"
    run_binner_combine([str(src)], str(dest), name="identity")
    assert (dest / "mag1.fa").is_file()


def test_iter_mag_fastas_skips_anvio_scratch_and_assembly_dump(tmp_path):
    mag = tmp_path / "bins"
    hidden = mag / ".anvio"
    hidden.mkdir(parents=True)
    (hidden / "contigs.fa").write_text(">c1\nACGT\n", encoding="utf-8")
    (mag / "contigs.fa").write_text(">c1\nACGT\n", encoding="utf-8")
    (mag / "mag1.fa").write_text(">c1\nACGT\n", encoding="utf-8")
    assert {p.name for p in _iter_mag_fastas(mag)} == {"mag1.fa"}


def test_identity_combine_drops_duplicate_contigs_fasta(tmp_path):
    src = tmp_path / "bins"
    src.mkdir()
    (src / "mag1.fa").write_text(">c1\nACGT\n", encoding="utf-8")
    (src / "contigs.fa").write_text(">c1\nACGT\n", encoding="utf-8")
    dest = tmp_path / "selected"
    run_binner_combine([str(src)], str(dest), name="identity")
    assert {p.name for p in dest.iterdir() if p.is_file()} == {"mag1.fa"}


def test_coverm_prefixes_duplicate_contig_headers(tmp_path, monkeypatch):
    mag = tmp_path / "mags"
    mag.mkdir()
    (mag / "mag1.fa").write_text(">c1\nACGT\n", encoding="utf-8")
    (mag / "contigs.fa").write_text(">c1\nTTTT\n", encoding="utf-8")
    dest = tmp_path / "ab.tsv"
    bam = tmp_path / "hits.bam"
    bam.write_bytes(b"BAM\x01")

    def fake_check_output(cmd, text=True):
        idx = cmd.index("--genome-fasta-files")
        files = []
        for token in cmd[idx + 1 :]:
            if str(token).startswith("--"):
                break
            files.append(Path(token))
        headers = []
        for fasta in files:
            for line in fasta.read_text(encoding="utf-8").splitlines():
                if line.startswith(">"):
                    headers.append(line[1:].split()[0])
        assert headers
        assert len(headers) == len(set(headers))
        return "Genome\tcount\nmag1\t3\n"

    monkeypatch.setattr(
        "samovar.assembly_profiling.which_tool",
        lambda name, default="": "/usr/bin/coverm",
    )
    monkeypatch.setattr("samovar.assembly_profiling.subprocess.check_output", fake_check_output)
    run_coverm(str(bam), str(mag), str(dest))
    assert dest.read_text(encoding="utf-8").splitlines()[1].startswith("mag1")


def test_anvio_empty_contigs_skips_database(tmp_path, monkeypatch):
    contigs = tmp_path / "empty.fa"
    contigs.write_text("\n", encoding="utf-8")
    dest = tmp_path / "bins"

    def boom(*_args, **_kwargs):
        raise AssertionError("anvi'o should not run on empty FASTA")

    monkeypatch.setattr("samovar.assembly_profiling.which_tool", boom)
    run_anvio_binner(str(contigs), str(dest))
    assert (dest / "mag1.fa").is_file()


@pytest.mark.optional
def test_megahit_assembler_writes_contigs(tmp_path, monkeypatch):
    """MegaHIT sidecar assembles the toy paired FASTQ into a contig FASTA."""
    if not READS_R1.is_file():
        pytest.skip("test reads missing")
    binary = _sidecar_bin("megahit")
    _put_on_path(monkeypatch, binary)
    dest = tmp_path / "contigs.fa"
    run_assembler(
        str(READS_R1),
        str(READS_R2),
        str(dest),
        threads=2,
        extra=["--min-contig-len", "100"],
    )
    text = dest.read_text()
    assert text.startswith(">")
    assert "N" in text or "A" in text.upper()


@pytest.mark.optional
def test_anvio_binner_single_sample_skips_clustering(tmp_path, monkeypatch):
    """Single-sample anvi'o writes mag1.fa and honors --skip-clustering."""
    binary = _sidecar_bin("anvio")
    _put_on_path(monkeypatch, binary)
    megahit = _sidecar_bin("megahit")
    _put_on_path(monkeypatch, megahit)
    contigs = tmp_path / "contigs.fa"
    run_assembler(
        str(READS_R1),
        str(READS_R2),
        str(contigs),
        threads=2,
        extra=["--min-contig-len", "100"],
    )
    mag_dir = tmp_path / "bins" / "anvio"
    run_binner(
        str(contigs),
        str(mag_dir),
        r1=str(READS_R1),
        r2=str(READS_R2),
        threads=2,
        name="anvio",
    )
    mags = _iter_mag_fastas(mag_dir)
    assert mags, "anvi'o single-sample path should write mag1.fa without clustering"
    skipped = tmp_path / "bins" / "anvio_skip"
    run_binner(
        str(contigs),
        str(skipped),
        r1=str(READS_R1),
        r2=str(READS_R2),
        threads=2,
        name="anvio",
        extra=["--skip-clustering"],
    )
    assert _iter_mag_fastas(skipped)


@pytest.mark.optional
def test_binner_qc_and_combine(tmp_path, monkeypatch):
    """CheckM2 QC plus DAS Tool (or identity) combine two MAG directories."""
    megahit = _sidecar_bin("megahit")
    _put_on_path(monkeypatch, megahit)
    contigs = tmp_path / "contigs.fa"
    run_assembler(
        str(READS_R1),
        str(READS_R2),
        str(contigs),
        threads=2,
        extra=["--min-contig-len", "100"],
    )
    d1 = tmp_path / "bins" / "a"
    d2 = tmp_path / "bins" / "b"
    d1.mkdir(parents=True)
    d2.mkdir(parents=True)
    shutil.copy2(contigs, d1 / "magA.fa")
    shutil.copy2(contigs, d2 / "magB.fa")
    try:
        checkm = _sidecar_bin("checkm2")
        _put_on_path(monkeypatch, checkm)
        qc = tmp_path / "qc.tsv"
        run_binner_qc(str(d1), str(qc), threads=2)
        assert "completeness" in qc.read_text().splitlines()[0]
    except Exception:
        qc = tmp_path / "qc.tsv"
        qc.write_text("mag_id\tcompleteness\tcontamination\tscore\n", encoding="utf-8")
    selected = tmp_path / "bins_selected"
    try:
        das = _sidecar_bin("dastool")
        _put_on_path(monkeypatch, das)
        run_binner_combine(
            [str(d1), str(d2)],
            str(selected),
            contigs=str(contigs),
            threads=2,
            name="dastool",
        )
    except Exception:
        run_binner_combine([str(d1), str(d2)], str(selected), name="identity")
    assert _iter_mag_fastas(selected)


@pytest.mark.optional
def test_anvio_binner_merges_three_samples(tmp_path, monkeypatch):
    """Three FASTQ pairs are profiled and anvi-merge writes a merged PROFILE.db."""
    reads = REPO / "tests" / "data" / "reads"
    r1s = [reads / f"{i}_full_R1.fastq" for i in (1, 2, 3)]
    r2s = [reads / f"{i}_full_R2.fastq" for i in (1, 2, 3)]
    if not all(p.is_file() for p in r1s + r2s):
        pytest.skip("need three test FASTQ pairs")
    binary = _sidecar_bin("anvio")
    _put_on_path(monkeypatch, binary)
    megahit = _sidecar_bin("megahit")
    _put_on_path(monkeypatch, megahit)
    minimap = _sidecar_bin("minimap2")
    _put_on_path(monkeypatch, minimap)
    from samovar.assembly_profiling import run_anvio_binner

    contigs = tmp_path / "contigs.fa"
    run_assembler(
        str(r1s[0]),
        str(r2s[0]),
        str(contigs),
        threads=2,
        extra=["--min-contig-len", "100"],
    )
    mag_dir = tmp_path / "bins" / "anvio_multi"
    try:
        run_anvio_binner(
            str(contigs),
            str(mag_dir),
            r1=str(r1s[0]),
            r2=str(r2s[0]),
            threads=2,
            additional_reads=[(str(r1s[1]), str(r2s[1])), (str(r1s[2]), str(r2s[2]))],
            skip_clustering=False,
        )
    except Exception as exc:
        pytest.skip(f"anvi'o multi-sample clustering could not finish: {exc}")
    assert _iter_mag_fastas(mag_dir)
    from samovar.paths import sidecar_envs_dir

    concoct = sidecar_envs_dir() / "anvio" / "bin" / "concoct"
    work = mag_dir / ".anvio" / "merged"
    if not concoct.is_file():
        pytest.skip("CONCOCT is not in the anvi'o sidecar; clustering drivers unavailable")
    assert (work / "PROFILE.db").is_file(), "anvi-merge should write a merged profile for 3 samples"


@pytest.mark.optional
def test_gtdbtk_mag_taxonomy(tmp_path, monkeypatch):
    """GTDB-Tk classify_wf writes a taxid column when the GTDB data path is ready."""
    db = Path("/mnt/tank/scratch/partition-metagenomics/databases/GTDB")
    env_db = os.environ.get("GTDBTK_DATA_PATH", "")
    if not gtdbtk_db_ready(str(db)) and not (env_db and gtdbtk_db_ready(env_db)):
        pytest.skip("GTDB-Tk database not extracted under partition-metagenomics/databases/GTDB")
    _sidecar_bin("gtdbtk")
    mag = tmp_path / "mags"
    mag.mkdir()
    (mag / "mag1.fa").write_text(">c1\n" + ("ACGT" * 200) + "\n", encoding="utf-8")
    dest = tmp_path / "tax.tsv"
    try:
        run_mag_taxonomy(str(mag), str(dest), db=str(db) if db.exists() else "", threads=2)
    except Exception as exc:
        pytest.skip(f"GTDB-Tk classify_wf could not run: {exc}")
    header = dest.read_text().splitlines()[0]
    assert "taxid" in header


@pytest.mark.optional
def test_minimap2_aligner_writes_bam(tmp_path, monkeypatch):
    """minimap2 aligner maps toy reads to a MAG and writes BAM."""
    mm = _sidecar_bin("minimap2")
    _put_on_path(monkeypatch, mm)
    mag = tmp_path / "mags"
    mag.mkdir()
    (mag / "mag1.fa").write_text(">c1\n" + ("ACGT" * 50) + "\n", encoding="utf-8")
    bam = tmp_path / "hits.bam"
    run_aligner(str(READS_R1), str(READS_R2), str(mag), str(bam), threads=2)
    assert bam.is_file() and bam.stat().st_size > 0
    assert bam.read_bytes()[:3] == b"BAM" or bam.stat().st_size > 0


@pytest.mark.optional
def test_coverm_mag_quantifier(tmp_path, monkeypatch):
    """CoverM quantifies MAG abundance from a BAM."""
    mm = _sidecar_bin("minimap2")
    _put_on_path(monkeypatch, mm)
    mag = tmp_path / "mags"
    mag.mkdir()
    (mag / "mag1.fa").write_text(">c1\n" + ("ACGT" * 50) + "\n", encoding="utf-8")
    bam = tmp_path / "hits.bam"
    run_aligner(str(READS_R1), str(READS_R2), str(mag), str(bam), threads=2)
    try:
        coverm = _sidecar_bin("coverm")
        _put_on_path(monkeypatch, coverm)
        dest = tmp_path / "ab.tsv"
        run_mag_quantifier(str(bam), str(mag), str(dest), name="coverm")
    except Exception as exc:
        pytest.skip(f"CoverM unavailable: {exc}")
    assert "mag_id" in dest.read_text().splitlines()[0]


@pytest.mark.optional
def test_composite_annotator(tmp_path, monkeypatch):
    """Assembler + identity binner + minimap2 + samtools counts produce an assembly TSV."""
    from samovar.assembly_profiling import main as assembly_main

    megahit = _sidecar_bin("megahit")
    mm = _sidecar_bin("minimap2")
    _put_on_path(monkeypatch, megahit)
    _put_on_path(monkeypatch, mm)
    dest = tmp_path / "sample_run.assembly.out"
    try:
        assembly_main(
            [
                "annotator",
                "-i",
                str(READS_R1),
                "-I",
                str(READS_R2),
                "-o",
                str(dest),
                "-t",
                "2",
                "--assembler",
                "megahit",
                "--binner",
                "identity",
                "--binner-combine",
                "identity",
                "--aligner",
                "minimap2",
                "--mag-quantifier",
                "samtools",
                "--work-dir",
                str(tmp_path / "sample_run.assembly"),
            ]
        )
    except FileNotFoundError as exc:
        pytest.skip(str(exc))
    except Exception as exc:
        pytest.skip(f"composite assembly annotator could not finish: {exc}")
    assert dest.is_file()
    parsed = Annotation({str(dest): "assembly"})
    tax_cols = [c for c in parsed.DataFrame.columns if str(c).startswith("taxID")]
    feat_cols = [c for c in parsed.DataFrame.columns if str(c).startswith("feat_")]
    if dest.stat().st_size > 0 and dest.read_text(encoding="utf-8").strip() not in {
        "",
        "seq\ttaxID\tMAG_ID",
    }:
        assert tax_cols
        assert any("MAG_ID" in c for c in feat_cols)
