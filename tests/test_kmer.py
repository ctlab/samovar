"""Native k-mer cases are optional; shared feature/scoring contracts are core."""

import pandas as pd
import pytest

from samovar.annotation_columns import (
    feat_annotator_columns,
    is_tax_column,
    prefix_tool_columns,
    select_scoring_annotators,
)
from samovar.annotators_wrapper import get_annotator_instance
from samovar.kmer2 import main, parse_output
from samovar.parse_annotators import Annotation
from samovar.reprofiling import preprocess_data, train_models
from samovar.scores import score_annotators



def _fastq(records):
    chunks = []
    for read_id, seq in records:
        chunks.append(f"@{read_id}\n{seq}\n+\n{'I' * len(seq)}\n")
    return "".join(chunks)


@pytest.mark.optional
def test_kmer_counts_dinucleotides(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text(_fastq([("r0", "AAAT")]))
    out = tmp_path / "kmer.out"
    assert main(["-i", str(r1), "-o", str(out), "-t", "1"]) == 0
    frame = parse_output(str(out)).set_index("seq")
    row = frame.loc["r0"]
    assert len(frame.columns) == 16
    assert int(row["AA"]) == 2
    assert int(row["AT"]) == 1
    assert int(row["TA"]) == 0
    assert sum(int(row[c]) for c in frame.columns) == 3


@pytest.mark.optional
def test_kmer_skips_non_acgt(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text(_fastq([("r0", "ANAA")]))
    out = tmp_path / "kmer.out"
    assert main(["-i", str(r1), "-o", str(out), "-t", "1"]) == 0
    row = parse_output(str(out)).set_index("seq").loc["r0"]
    assert int(row["AA"]) == 1
    assert sum(int(row[c]) for c in row.index) == 1


@pytest.mark.optional
def test_kmer_cli_and_parse_output(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    r1.write_text(_fastq([("r0", "AAAT")]))
    r2.write_text(_fastq([("r0", "GGGG")]))
    out = tmp_path / "kmer.out"
    assert main(["-i", str(r1), "-I", str(r2), "-o", str(out), "-t", "1"]) == 0
    text = out.read_text()
    assert text.startswith("seq\tAA\t")
    frame = parse_output(str(out))
    assert "seq" in frame.columns
    assert "AA" in frame.columns
    assert "taxID" not in frame.columns
    row = frame.set_index("seq").loc["r0"]
    assert int(row["AA"]) == 2
    assert int(row["GG"]) == 3
    labeled = prefix_tool_columns(frame, "kmer2", 0)
    feat = feat_annotator_columns(labeled.columns)
    assert len(feat) == 16
    assert "feat_kmer2_0_AA" in labeled.columns
    assert not any(is_tax_column(c) for c in labeled.columns)


@pytest.mark.optional
def test_kmer_factory_and_annotation_class(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text("@readA|taxid:562|\nATGC\n+\nIIII\n")
    out = tmp_path / "sample_kmer2.kmer2.out"
    inst = get_annotator_instance("kmer2", {"run_name": "kmer2", "cmd": "kmer2"}, {})
    assert inst.__class__.__name__ == "Kmer2Annotator"
    assert main(["-i", str(r1), "-I", "", "-o", str(out), "-t", "1"]) == 0
    kaiju = tmp_path / "sample_kaiju.kaiju.out"
    kaiju.write_text("C\treadA|taxid:562|\t562\n")
    ann = Annotation({str(kaiju): "kaiju", str(out): "kmer2"})
    assert list(ann.tr().columns)
    assert len(ann.ft().columns) == 16
    assert all(c.startswith("feat_") for c in ann.ft().columns)


@pytest.mark.optional
def test_kmer_max_keeps_highest_variance(tmp_path):
    """Top k-mers by across-read variance, capped by the real k-mer count."""
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text(
        _fastq(
            [
                ("r1", "AAAA"),
                ("r2", "AAAA"),
                ("r3", "ACAC"),
                ("r4", "ACAC"),
            ]
        )
    )
    out = tmp_path / "kmer.out"
    assert main(["-i", str(r1), "-o", str(out), "-t", "2", "--kmer-max", "2"]) == 0
    frame = parse_output(str(out))
    assert list(frame.columns) == ["seq", "AA", "AC"]
    assert main(["-i", str(r1), "-o", str(out), "-t", "1", "--kmer-max", "100"]) == 0
    wide = parse_output(str(out))
    assert len([c for c in wide.columns if c != "seq"]) == 16


@pytest.mark.optional
def test_kmer_threads_and_duplicates_match(tmp_path):
    records = [("r%03d" % i, "ACGT" * (3 + (i % 5)) + ("N" if i % 7 == 0 else "")) for i in range(40)]
    records.append(("r000", "TTTT"))
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    r1.write_text(_fastq(records))
    r2.write_text(_fastq([("r000", "aaat"), ("only2", "GG")]))
    one = tmp_path / "t1.out"
    many = tmp_path / "t8.out"
    assert main(["-i", str(r1), "-I", str(r2), "-o", str(one), "-t", "1"]) == 0
    assert main(["-i", str(r1), "-I", str(r2), "-o", str(many), "-t", "8", "--kmer-max", "16"]) == 0
    assert one.read_text() == many.read_text()
    frame = parse_output(str(one)).set_index("seq")
    assert int(frame.loc["r000", "AA"]) == 2
    assert int(frame.loc["r000", "AT"]) == 1
    assert int(frame.loc["r000", "TT"]) == 3
    assert int(frame.loc["only2", "GG"]) == 1
    assert "only2" in frame.index


def test_scoring_ignores_feat_columns():
    df = pd.DataFrame(
        {
            "taxID_kaiju_0": ["562", "562", "9606"],
            "feat_kmer2_0_AA": [9, 8, 1],
            "true": ["562", "562", "9606"],
        }
    )
    table = score_annotators(df, ["taxID_kaiju_0", "feat_kmer2_0_AA"])
    names = set(table["annotator"].astype(str).str.lower())
    assert not any("feat" in n for n in names)
    assert select_scoring_annotators(df, ["feat_kmer2_0_AA", "taxID_kaiju_0"]) == [
        "taxID_kaiju_0"
    ]
    mag_df = pd.DataFrame(
        {
            "taxID_assembly_0": ["562", "562"],
            "feat_assembly_0_MAG_ID": ["mag1", "mag1"],
            "true": ["562", "562"],
        }
    )
    assert select_scoring_annotators(mag_df) == ["taxID_assembly_0"]


def test_ml_uses_feat_columns():
    df = pd.DataFrame(
        {
            "seq": [f"s{i}" for i in range(16)],
            "taxID_kaiju_0": ["1", "2"] * 8,
            "feat_kmer2_0_AA": list(range(16)),
            "length": [4] * 16,
            "true": [1, 2] * 8,
        }
    )
    processed = preprocess_data(df)
    assert "feat_kmer2_0_AA" in processed.columns
    best, models, metrics, feature_cols = train_models(processed, test_size=0.25)
    assert "feat_kmer2_0_AA" in feature_cols
    assert any(c.startswith("taxid_") for c in feature_cols)
    assert best is not None
