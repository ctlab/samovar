import gzip
import bz2
import json
import os
import pickle
from pathlib import Path
from unittest.mock import patch
import pytest

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from samovar.build_database import (
    build_database_from_config,
    find_local_proteome_or_gff,
    kaiju_protein_id,
    nucleotide_to_orf_records,
    process_fasta_kaiju,
    process_fasta_kraken2,
    translate_cds_from_gff,
)


def test_metaphlan_build_checks_marker_pair_and_indexes_via_shared_runner(tmp_path, monkeypatch):
    import yaml
    fasta = tmp_path / "markers.fna"
    fasta.write_text(">SGB1__gene\nACGT\n")
    info = tmp_path / "markers.pkl"
    with bz2.open(info, "wb") as handle:
        pickle.dump({"taxonomy": {"s__example": ("562", 4)},
                     "markers": {"SGB1__gene": {"len": 4}}}, handle)
    cfg = tmp_path / "build.yaml"
    cfg.write_text(yaml.safe_dump(dict(marker_fasta=str(fasta), marker_metadata=str(info))))
    dest = tmp_path / "db with spaces"
    calls = []

    def build(cmd):
        calls.append(cmd)
        assert cmd[:3] == ["bowtie2-build", "--threads", "2"]
        assert Path(cmd[-2]).read_text() == fasta.read_text()
        for part in ("1", "2", "3", "4", "rev.1", "rev.2"):
            Path(f"{cmd[-1]}.{part}.bt2l").write_bytes(b"index")

    monkeypatch.setattr("samovar.build_database.run_command", build)
    build_database_from_config(str(cfg), db_type="metaphlan", db_path=str(dest), threads=2)
    manifest = json.loads((dest / "metaphlan.manifest.json").read_text())
    assert manifest["marker_count"] == 1 and len(manifest["fasta_sha256"]) == 64
    assert (dest / "markers.pkl").read_bytes() == info.read_bytes()
    assert len(calls) == 1
    fasta.write_text(">SGB1__gene\nAC\n")
    with pytest.raises(ValueError, match="length"):
        build_database_from_config(str(cfg), db_type="metaphlan", db_path=str(dest), threads=2)
    assert len(calls) == 1



def test_sylph_build_preserves_ncbi_mapping_and_resources(tmp_path, monkeypatch):
    import yaml
    genomes = tmp_path / "genomes with spaces"
    genomes.mkdir()
    (genomes / "562.fna").write_text(">contig\nACGTACGT\n")
    cfg = tmp_path / "build.yaml"
    cfg.write_text(yaml.safe_dump(dict(input_dir=[str(genomes)], compression=50)))
    dest = tmp_path / "db with spaces"
    calls = []

    def sketch(cmd):
        calls.append(cmd)
        assert cmd[:3] == ["sylph", "sketch", "-l"]
        assert Path(cmd[3]).read_text().splitlines() == [str(genomes / "562.fna")]
        assert cmd[-4:] == ["-t", "3", "-c", "50"]
        Path(cmd[cmd.index("-o") + 1] + ".syldb").write_bytes(b"index")

    monkeypatch.setattr("samovar.build_database.run_command", sketch)
    build_database_from_config(str(cfg), db_type="sylph", db_path=str(dest), threads=3)
    assert (dest / "genomes.syldb").read_bytes() == b"index"
    assert (dest / "genomes.syldb.taxonomy.tsv").read_text().splitlines() == [
        "Genome_file\ttaxid", str(genomes / "562.fna") + "\t562"]
    manifest = json.loads((dest / "sylph.manifest.json").read_text())
    assert manifest["compression"] == 50 and manifest["genomes"][0]["taxid"] == "562"
    assert len(calls) == 1


def test_sylph_build_rejects_missing_index_and_unresolved_taxonomy(tmp_path, monkeypatch):
    from samovar.build_database import build_database_sylph
    genome = tmp_path / "genome.fna"
    genome.write_text(">contig\nACGTACGT\n")
    dest = tmp_path / "db"
    monkeypatch.setattr("samovar.build_database.run_command", lambda cmd: None)
    with pytest.raises(ValueError, match="positive NCBI taxid"):
        build_database_sylph({str(genome): "SGB1"}, str(dest))
    with pytest.raises(RuntimeError, match="did not create"):
        build_database_sylph({str(genome): "562"}, str(dest))
    assert not (dest / "genomes.syldb.taxonomy.tsv").exists()
    genome.write_text(">protein\nMKKLL\n")
    with pytest.raises(ValueError, match="No nucleotide genomes"):
        build_database_sylph({str(genome): "562"}, str(dest))

def test_kaiju_protein_id_puts_taxid_last():
    assert kaiju_protein_id("9606") == "9606"
    assert kaiju_protein_id("9606", "NP_1").endswith("_9606")
    assert kaiju_protein_id("562", "acc").rsplit("_", 1)[-1] == "562"


def _write_fasta(path, records):
    SeqIO.write(records, path, "fasta")


def _taxid_from_kaiju_header(seq_id: str) -> str:
    return seq_id.rsplit("_", 1)[-1]


def test_process_fasta_kaiju_protein_headers_end_with_taxid(tmp_path):
    faa = tmp_path / "562.faa"
    _write_fasta(
        faa,
        [SeqRecord(Seq("MEEITAKY"), id="NP_1", description="ecoli protein")],
    )
    out = process_fasta_kaiju(str(faa), "562", protein=True, fetch_missing=False, min_aa=3)
    recs = list(SeqIO.parse(out, "fasta"))
    os.remove(out)
    assert recs
    assert _taxid_from_kaiju_header(recs[0].id) == "562"
    assert "*" not in str(recs[0].seq)


def test_process_fasta_kaiju_orfs_from_nucleotide(tmp_path):
    # ATG GAA TTC GGT TAA  -> MEFG (plus other-frame ORFs)
    fna = tmp_path / "562.fna"
    _write_fasta(fna, [SeqRecord(Seq("ATGGAATTCGGTTAA"), id="gene1")])
    out = process_fasta_kaiju(str(fna), "562", fetch_missing=False, min_aa=3)
    recs = list(SeqIO.parse(out, "fasta"))
    os.remove(out)
    assert recs
    assert all(_taxid_from_kaiju_header(rec.id) == "562" for rec in recs)
    assert all("*" not in str(rec.seq) for rec in recs)
    assert any("MEFG" in str(rec.seq) or str(rec.seq) == "MEFG" for rec in recs)


def test_translate_cds_from_gff(tmp_path):
    fna = tmp_path / "9606.fna"
    gff = tmp_path / "9606.gff"
    _write_fasta(fna, [SeqRecord(Seq("ATGGAATTCTAA"), id="chr1")])
    gff.write_text(
        "chr1\tsrc\tCDS\t1\t12\t.\t+\t0\tID=cds1;protein_id=prot1\n"
    )
    recs = translate_cds_from_gff(str(fna), str(gff), "9606", min_aa=2)
    assert len(recs) == 1
    assert _taxid_from_kaiju_header(recs[0].id) == "9606"
    assert str(recs[0].seq) == "MEF"


def test_process_fasta_kaiju_uses_sibling_faa(tmp_path):
    fna = tmp_path / "123.fna"
    faa = tmp_path / "123.faa"
    _write_fasta(fna, [SeqRecord(Seq("AAAAAA"), id="genome")])
    _write_fasta(faa, [SeqRecord(Seq("MKKLL"), id="prot")])
    kind, path = find_local_proteome_or_gff(str(fna), "123")
    assert kind == "protein"
    assert path == str(faa)
    out = process_fasta_kaiju(str(fna), "123", fetch_missing=False, min_aa=3)
    recs = list(SeqIO.parse(out, "fasta"))
    os.remove(out)
    assert _taxid_from_kaiju_header(recs[0].id) == "123"
    assert str(recs[0].seq) == "MKKLL"


def test_process_fasta_kaiju_does_not_fetch_when_disabled(tmp_path):
    fna = tmp_path / "999.fna"
    _write_fasta(fna, [SeqRecord(Seq("ATGAAATAA"), id="x")])
    with patch("samovar.build_database.fetch_proteome") as fetch_p, patch(
        "samovar.build_database.fetch_gff"
    ) as fetch_g:
        process_fasta_kaiju(str(fna), "999", fetch_missing=False, min_aa=2)
        fetch_p.assert_not_called()
        fetch_g.assert_not_called()


def test_nucleotide_to_orf_records_min_length(tmp_path):
    fna = tmp_path / "1.fna"
    _write_fasta(fna, [SeqRecord(Seq("ATGAAATAA"), id="short")])  # MK*
    recs = nucleotide_to_orf_records(str(fna), "1", min_aa=10)
    assert recs == []


def test_process_fasta_kraken2_gzipped_input(tmp_path):
    fna = tmp_path / "562.fna.gz"
    with gzip.open(fna, "wt") as handle:
        SeqIO.write([SeqRecord(Seq("ACGTACGT"), id="c1")], handle, "fasta")
    out = process_fasta_kraken2(str(fna), "562")
    recs = list(SeqIO.parse(out, "fasta"))
    os.remove(out)
    assert recs[0].id == "seq0|kraken:taxid|562"


def test_find_local_proteome_gzipped_faa(tmp_path):
    fna = tmp_path / "123.fna"
    faa = tmp_path / "123.faa.gz"
    _write_fasta(fna, [SeqRecord(Seq("AAAAAA"), id="genome")])
    with gzip.open(faa, "wt") as handle:
        SeqIO.write([SeqRecord(Seq("MKKLL"), id="prot")], handle, "fasta")
    kind, path = find_local_proteome_or_gff(str(fna), "123")
    assert kind == "protein"
    assert path == str(faa)


def test_process_fasta_kraken2_keeps_all_contigs(tmp_path):
    fna = tmp_path / "562.fna"
    _write_fasta(
        fna,
        [
            SeqRecord(Seq("ACGTACGT"), id="c1"),
            SeqRecord(Seq("GGGGCCCC"), id="c2"),
        ],
    )
    out = process_fasta_kraken2(str(fna), "562")
    recs = list(SeqIO.parse(out, "fasta"))
    os.remove(out)
    assert len(recs) == 2
    assert recs[0].id == "seq0|kraken:taxid|562"
    assert recs[1].id == "seq1|kraken:taxid|562"
