"""Builtin k-mer encoder Feature extractor (metauto encoding stage)."""

from pathlib import Path

import pytest

from samovar.annotation_columns import feat_annotator_columns, prefix_tool_columns
from samovar.annotators_wrapper import get_annotator_instance
from samovar.kmer_encoder import kmer_ids, main, parse_output, seq_to_vec
from samovar.parse_annotators import Annotation, match_annotation


def test_seq_to_vec_normalizes():
    vec = seq_to_vec("AAAT", k=2)
    ids = kmer_ids(2)
    by = dict(zip(ids, vec))
    assert pytest.approx(by["AA"] + by["AT"], 1e-9) == 1.0
    assert by["AA"] == pytest.approx(2 / 3)
    assert by["GG"] == 0


def test_kmer_encoder_freq_cli(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text("@r0\nAAAT\n+\nIIII\n")
    out = tmp_path / "enc.out"
    assert main(["-i", str(r1), "-I", "", "-o", str(out), "-k", "2"]) == 0
    text = out.read_text()
    assert text.startswith("seq\tAA\t")
    frame = parse_output(str(out))
    assert "seq" in frame.columns
    assert "taxID" not in frame.columns
    assert "AA" in frame.columns
    labeled = prefix_tool_columns(frame, "kmer_encoder", 0)
    feat = feat_annotator_columns(labeled.columns)
    assert len(feat) == 16
    assert any("AA" in c for c in feat)


def test_factory_and_match_annotation(tmp_path):
    inst = get_annotator_instance("kmer_encoder", {"run_name": "kmer-encoder-test"}, {})
    assert inst.__class__.__name__ == "KmerEncoderAnnotator"
    outs = inst.get_expected_outputs("1_full", "/tmp")
    assert outs[0].endswith("1_full_kmer-encoder-test.kmer_encoder.out")
    assert match_annotation("1_full_kmer-encoder-test.kmer_encoder.out") == "kmer_encoder"
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text("@readA|taxid:562|\nATGC\n+\nIIII\n")
    dest = tmp_path / "sample_kmer-encoder-test.kmer_encoder.out"
    assert main(["-i", str(r1), "-o", str(dest), "-k", "2"]) == 0
    kaiju = tmp_path / "sample_kaiju.kaiju.out"
    kaiju.write_text("C\treadA|taxid:562|\t562\n")
    ann = Annotation({str(kaiju): "kaiju", str(dest): "kmer_encoder"})
    assert list(ann.tr().columns)
    assert len(ann.ft().columns) == 16
    assert all(str(c).startswith("feat_") for c in ann.ft().columns)


def test_train_and_latent_features(tmp_path):
    pytest.importorskip("torch")
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "562.fna").write_text(">g\n" + ("ATGC" * 20) + "\n")
    (genomes / "9606.fna").write_text(">h\n" + ("GGCC" * 20) + "\n")
    model = tmp_path / "kmer_encoder.pt"
    assert main(
        ["train", "-c", str(genomes), "-o", str(model), "-k", "2", "--epochs", "1", "--latent", "8"]
    ) == 0
    assert model.is_file()
    r1 = tmp_path / "r1.fastq"
    r1.write_text("@r0\nATGCATGC\n+\nIIIIIIII\n")
    out = tmp_path / "z.out"
    assert main(["-i", str(r1), "-d", str(model), "-o", str(out)]) == 0
    frame = parse_output(str(out))
    zcols = [c for c in frame.columns if str(c).startswith("z")]
    assert len(zcols) == 8
    assert "taxID" not in frame.columns
    labeled = prefix_tool_columns(frame, "kmer_encoder", 0)
    assert all(str(c).startswith("feat_") or c == "seq" for c in labeled.columns)
