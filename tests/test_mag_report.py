"""MAG loss and MAG_ID shares in the MultiQC plot folders."""
import json

from samovar.MAG.report import mag_id_shares, subsample_lines, write_mag_mapping_plots

import pytest

# Detailed MAG report layout.
pytestmark = pytest.mark.optional


def test_large_mag_table_is_subsampled(tmp_path):
    path = tmp_path / "reads.tsv"
    lines = ["seq\ttaxID\tMAG_ID\n"]
    lines.extend(f"r{i}\t562\t{'M1' if i % 2 == 0 else 'M2'}\n" for i in range(400))
    path.write_text("".join(lines), encoding="utf-8")
    header, picked = subsample_lines(path, max_rows=20, full_scan_bytes=10)
    assert header.startswith("seq")
    assert len(picked) == 20
    shares = mag_id_shares(path, max_rows=20, full_scan_bytes=10)
    assert shares[0][0] == ""
    assert set(shares[0][1]) <= {"M1", "M2"}
    assert abs(sum(shares[0][1].values()) - 100) < 0.1


def test_combined_and_reprofiled_tables_keep_mag_sections(tmp_path):
    regen = tmp_path / "regenerated_annotations"
    regen.mkdir()
    (regen / "combined_annotation_table.csv").write_text(
        "seq,sample,feat_assembly_0_MAG_ID\n"
        "a,s1,M1\n"
        "b,s1,\n"
        "c,s2,M2\n"
        "d,s2,M2\n",
        encoding="utf-8",
    )
    repro = tmp_path / "reprofiled_annotations"
    repro.mkdir()
    (repro / "s1_reprofiled.csv").write_text(
        "seq,feat_assembly_0_MAG_ID\na,M1\nb,\n",
        encoding="utf-8",
    )
    written = write_mag_mapping_plots(tmp_path)
    names = {path.name for path in written}
    assert "mag_regenerated_ids_mqc.json" in names
    assert "mag_reprofiled_ids_mqc.json" in names
    assert "mag_reprofiled_losses_mqc.json" in names
    regen_ids = json.loads(
        (tmp_path / "regenerated_annotations_plots" / "mag_regenerated_ids_mqc.json").read_text()
    )
    assert regen_ids["data"]["s1 feat_assembly_0"]["M1"] == 50.0
    assert regen_ids["data"]["s1 feat_assembly_0"]["(unassigned)"] == 50.0
    assert regen_ids["data"]["s2 feat_assembly_0"]["M2"] == 100.0
    repro_ids = json.loads(
        (tmp_path / "reprofiled_annotations_plots" / "mag_reprofiled_ids_mqc.json").read_text()
    )
    assert repro_ids["data"]["s1 feat_assembly_0"]["(unassigned)"] == 50.0
    losses = json.loads(
        (tmp_path / "reprofiled_annotations_plots" / "mag_reprofiled_losses_mqc.json").read_text()
    )
    assert losses["data"]["s1 feat_assembly_0"]["unassigned"] == 50.0
    from samovar.stage_report import bundle_multiqc

    staged = bundle_multiqc(tmp_path)
    bundled = "\n".join(path.read_text(encoding="utf-8") for path in staged.glob("*mag_*"))
    assert "MAG_ID read shares" in bundled
    assert "unassigned" in bundled


def test_mapping_losses_and_reprofiled_mag_shares(tmp_path):
    reports = tmp_path / "initial_reports"
    reports.mkdir()
    (reports / "s1_assembly.assembly.out").write_text(
        "seq\ttaxID\tMAG_ID\n"
        "a\t562\tM1\n"
        "b\t562\tM1\n"
        "c\t0\t\n"
        "d\t9606\tM2\n",
        encoding="utf-8",
    )
    (reports / "s1_assembly.assembly.out.mapping.json").write_text(
        json.dumps(
            {
                "total": 4,
                "unmapped": 1,
                "multi_hit": 1,
                "unresolved": 0,
                "unclassified": 0,
                "unmapped_fraction": 0.25,
                "multi_hit_fraction": 0.25,
                "unresolved_fraction": 0,
                "unclassified_fraction": 0,
            }
        ),
        encoding="utf-8",
    )
    (tmp_path / "initial_annotations").mkdir()
    (tmp_path / "initial_annotations" / "s1.annotation.csv").write_text(
        "seq,taxID\na,562\n", encoding="utf-8"
    )
    repro = tmp_path / "reprofiled_annotations"
    repro.mkdir()
    (repro / "s1.annotation.csv").write_text(
        "seq,feat_assembly_0_MAG_ID\na,M1\nb,M1\nc,\nd,M2\n",
        encoding="utf-8",
    )
    written = write_mag_mapping_plots(tmp_path)
    names = {path.name for path in written}
    assert "mag_initial_losses_mqc.json" in names
    assert "mag_initial_ids_mqc.json" in names
    assert "mag_reprofiled_ids_mqc.json" in names
    losses = json.loads((tmp_path / "initial_annotations_plots" / "mag_initial_losses_mqc.json").read_text())
    assert losses["plot_type"] == "bargraph"
    assert losses["data"]["s1 assembly"]["unmapped"] == 25.0
    assert losses["data"]["s1 assembly"]["multi_hit"] == 25.0
    assert losses["data"]["s1 assembly"]["unresolved"] == 0.0
    ids = json.loads((tmp_path / "initial_annotations_plots" / "mag_initial_ids_mqc.json").read_text())
    assert ids["data"]["s1 assembly"]["M1"] == 50.0
    assert ids["data"]["s1 assembly"]["(unassigned)"] == 25.0
    repro_ids = json.loads((tmp_path / "reprofiled_annotations_plots" / "mag_reprofiled_ids_mqc.json").read_text())
    assert repro_ids["data"]["s1 feat_assembly_0"]["M1"] == 50.0
    assert repro_ids["data"]["s1 feat_assembly_0"]["(unassigned)"] == 25.0
    from samovar.stage_report import bundle_multiqc

    staged = bundle_multiqc(tmp_path)
    bundled = "\n".join(path.read_text(encoding="utf-8") for path in staged.glob("*mag_initial_losses*"))
    assert "multi_hit" in bundled
    assert "unresolved" in bundled
