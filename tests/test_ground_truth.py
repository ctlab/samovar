"""CAMI-style and from-genomes ground-truth table helpers."""

from pathlib import Path

import pytest

from samovar.ground_truth import (
    GROUND_TRUTH_TABLE,
    PARSE_GENOME,
    append_from_fasta,
    append_from_genome_dir,
    iter_truth_pairs,
    main as gt_main,
    normalize_regenerated_mode,
    taxid_from_header,
    write_normalized_table,
)

REPO = Path(__file__).resolve().parents[1]
CAMI_TOY = REPO / "tests" / "data" / "cami_toy"
TEST_GENOMES = REPO / "data" / "test_genomes"


def test_normalize_regenerated_mode_aliases():
    assert normalize_regenerated_mode(None) == PARSE_GENOME
    assert normalize_regenerated_mode("parse_genome") == PARSE_GENOME
    assert normalize_regenerated_mode("cami") == GROUND_TRUTH_TABLE
    assert normalize_regenerated_mode("ground-truth-table") == GROUND_TRUTH_TABLE
    with pytest.raises(ValueError):
        normalize_regenerated_mode("headers-please")


def test_iter_truth_pairs_cami_and_simple(tmp_path):
    cami = tmp_path / "map.tsv"
    cami.write_text(
        "#anonymous_read_id\tgenome_id\ttax_id\tread_id\n"
        "S0R0/1\tASV153.8\t818\tNZ_JBDQBI010000009.1-9644/1\n"
        "S0R0/2\tASV153.8\t818\tNZ_JBDQBI010000009.1-9644/2\n"
    )
    pairs = list(iter_truth_pairs(cami))
    keys = {k for k, _ in pairs}
    assert "S0R0/1" in keys
    assert ("S0R0/1", "818") in pairs
    simple = tmp_path / "simple.tsv"
    simple.write_text("seq\ttaxid\ncontigA\t562\n")
    assert list(iter_truth_pairs(simple)) == [("contigA", "562")]


def test_append_from_fasta_uses_filename_taxid(tmp_path):
    fa = tmp_path / "562.fna"
    fa.write_text(">NC_000913.3 Escherichia coli\nATGC\n>plasmid extra\nGGCC\n")
    dest = tmp_path / "gt.tsv"
    n = append_from_fasta(fa, dest)
    assert n == 2
    text = dest.read_text()
    assert text.splitlines()[0] == "seq\ttaxid"
    assert "NC_000913.3\t562" in text
    assert "plasmid\t562" in text


def test_from_genomes_cli_appends_test_genomes(tmp_path):
    dest = tmp_path / "gt.tsv"
    rc = gt_main(["from-genomes", "--genome-dir", str(TEST_GENOMES / "meta"), "-o", str(dest)])
    assert rc == 0
    rows = dest.read_text().strip().splitlines()
    assert rows[0] == "seq\ttaxid"
    assert any("\t562" in line for line in rows[1:])
    assert any("\t4932" in line for line in rows[1:])
    assert taxid_from_header("Ecoli.fna|taxid:562|-NC_000913.3", "0") == "562"


def test_from_genomes_overwrite(tmp_path):
    dest = tmp_path / "gt.tsv"
    dest.write_text("seq\ttaxid\nold\t1\n")
    fa_dir = tmp_path / "genomes"
    fa_dir.mkdir()
    (fa_dir / "9606.fna").write_text(">chr1 Homo\nAAAA\n")
    rc = gt_main(
        ["from-genomes", "--genome-dir", str(fa_dir), "-o", str(dest), "--overwrite"]
    )
    assert rc == 0
    text = dest.read_text()
    assert "old\t1" not in text
    assert "chr1\t9606" in text


def test_normalize_cami_toy_mapping_if_present():
    mapping = CAMI_TOY / "reads_mapping.tsv.gz"
    if not mapping.is_file():
        mapping = CAMI_TOY / "reads_mapping.tsv"
    if not mapping.is_file():
        pytest.skip("CAMI toy mapping not extracted")
    pairs = list(iter_truth_pairs(mapping))
    assert len(pairs) >= 1000
    tax = {t for _, t in pairs}
    assert tax
    assert "0" not in tax


def test_write_normalized_table(tmp_path):
    src = tmp_path / "cami.tsv"
    src.write_text("#anonymous_read_id\tgenome_id\ttax_id\tread_id\nA/1\tg\t818\tr/1\n")
    dest = tmp_path / "norm.tsv"
    assert write_normalized_table(src, dest) == 2
    assert dest.read_text().splitlines() == ["seq\ttaxid", "A/1\t818", "r/1\t818"]
