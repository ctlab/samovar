"""Builtin GC feature extractor."""

from samovar.annotation_columns import feat_annotator_columns, is_tax_column, prefix_tool_columns
from samovar.annotators_wrapper import get_annotator_instance
from samovar.gc import main, parse_output
from samovar.parse_annotators import Annotation, match_annotation
from samovar.qc import gc_fraction


def test_gc_fraction_matches_qc():
    assert abs(gc_fraction("ATGC") - 0.5) < 1e-9
    assert abs(gc_fraction("GGCC") - 1.0) < 1e-9


def test_gc_cli_and_prefix(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text("@r0\nATGC\n+\nIIII\n")
    out = tmp_path / "sample_gc.gc.out"
    assert main(["-i", str(r1), "-I", "", "-o", str(out), "-t", "1"]) == 0
    frame = parse_output(str(out))
    assert list(frame.columns)[:2] == ["seq", "GC"]
    assert "taxID" not in frame.columns
    labeled = prefix_tool_columns(frame, "gc", 0)
    feat = feat_annotator_columns(labeled.columns)
    assert feat == ["feat_gc_0_GC"]
    assert not any(is_tax_column(c) for c in labeled.columns)


def test_gc_factory_and_annotation(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text("@readA|taxid:562|\nGGGG\n+\nIIII\n")
    out = tmp_path / "sample_gc.gc.out"
    inst = get_annotator_instance("gc", {"run_name": "gc", "cmd": "gc"}, {})
    assert inst.__class__.__name__ == "GcAnnotator"
    assert main(["-i", str(r1), "-I", "", "-o", str(out)]) == 0
    kaiju = tmp_path / "sample_kaiju.kaiju.out"
    kaiju.write_text("C\treadA|taxid:562|\t562\n")
    ann = Annotation({str(kaiju): "kaiju", str(out): "gc"})
    assert list(ann.ft().columns) == ["feat_gc_1_GC"]
    assert match_annotation("s.gc.out") == "gc"
