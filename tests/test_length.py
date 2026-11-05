"""Builtin read-length feature extractor."""

from samovar.annotation_columns import feat_annotator_columns, is_tax_column, prefix_tool_columns
from samovar.annotators_wrapper import get_annotator_instance
from samovar.length import main, parse_output
from samovar.parse_annotators import Annotation, match_annotation


def test_length_cli_and_prefix(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    r1.write_text("@r0\nATGC\n+\nIIII\n")
    r2.write_text("@r0\nAAAAAA\n+\nIIIIII\n")
    out = tmp_path / "sample_length.length.out"
    assert main(["-i", str(r1), "-I", str(r2), "-o", str(out), "-t", "1"]) == 0
    frame = parse_output(str(out))
    assert list(frame.columns)[:2] == ["seq", "length"]
    assert int(frame.loc[0, "length"]) == 10
    assert "taxID" not in frame.columns
    labeled = prefix_tool_columns(frame, "length", 0)
    feat = feat_annotator_columns(labeled.columns)
    assert feat == ["feat_length_0_length"]
    assert not any(is_tax_column(c) for c in labeled.columns)


def test_length_factory_and_annotation(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r1.write_text("@readA|taxid:562|\nGGGG\n+\nIIII\n")
    out = tmp_path / "sample_length.length.out"
    inst = get_annotator_instance("length", {"run_name": "length", "cmd": "length"}, {})
    assert inst.__class__.__name__ == "LengthAnnotator"
    assert main(["-i", str(r1), "-I", "", "-o", str(out)]) == 0
    kaiju = tmp_path / "sample_kaiju.kaiju.out"
    kaiju.write_text("C\treadA|taxid:562|\t562\n")
    ann = Annotation({str(kaiju): "kaiju", str(out): "length"})
    assert list(ann.ft().columns) == ["feat_length_1_length"]
    assert match_annotation("s.length.out") == "length"
