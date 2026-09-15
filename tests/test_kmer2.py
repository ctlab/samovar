"""Builtin k=2 k-mer feature extractor and Feature contract helpers."""

from pathlib import Path

import pandas as pd

from samovar.annotation_columns import (
    feat_annotator_columns,
    is_feat_column,
    is_tax_column,
    prefix_tool_columns,
    select_scoring_annotators,
)
from samovar.annotators_wrapper import get_annotator_instance
from samovar.kmer2 import KMER_IDS, count_kmers, main, parse_output
from samovar.main_config import normalize_tool_group
from samovar.parse_annotators import Annotation
from samovar.reprofiling import preprocess_data, train_models
from samovar.scores import score_annotators


def test_feature_type_aliases_annotator():
    assert normalize_tool_group("feature") == "annotator"
    assert normalize_tool_group("seq2feature") == "annotator"


def test_kmer2_counts_dinucleotides():
    counts = count_kmers("AAAT")
    by = dict(zip(KMER_IDS, counts))
    assert len(KMER_IDS) == 16
    assert by["AA"] == 2
    assert by["AT"] == 1
    assert by["TA"] == 0
    assert sum(counts) == 3


def test_kmer2_cli_and_parse_output(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    rec = "@r0\nAAAT\n+\nIIII\n"
    r1.write_text(rec)
    r2.write_text("@r0\nGGGG\n+\nIIII\n")
    out = tmp_path / "kmer2.out"
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


def test_kmer2_factory_and_annotation_class(tmp_path):
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


def test_count_kmers_skips_non_acgt():
    # N-containing windows are skipped; AA from ANAA still counts once.
    assert count_kmers("ANAA")[KMER_IDS.index("AA")] == 1


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
