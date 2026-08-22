import gzip
import os
from pathlib import Path
from unittest.mock import patch

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from samovar.build_database import (
    find_local_proteome_or_gff,
    kaiju_protein_id,
    nucleotide_to_orf_records,
    process_fasta_kaiju,
    process_fasta_kraken2,
    translate_cds_from_gff,
)


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
