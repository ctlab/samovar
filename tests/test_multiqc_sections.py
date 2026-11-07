"""MultiQC section selection is rebuilt from a finished run."""

import json
from pathlib import Path

import pytest
import yaml

from samovar.multiqc_report import effective_settings, parse_cli, resolve_settings
from samovar.stage_report import bundle_multiqc


def _scores(root: Path) -> None:
    plots = root / "initial_annotations_plots"
    plots.mkdir(parents=True)
    (plots / "quality_scores.csv").write_text(
        "annotator,accuracy,f1\nkaiju,0.5,0.4\nkraken2,0.6,0.5\n",
        encoding="utf-8",
    )
    (plots / "F1_kaiju_mqc.json").write_text(
        json.dumps(
            {
                "plot_type": "heatmap",
                "section_name": "F1 heatmap — kaiju",
                "data": [[1]],
                "xcats": ["562"],
                "ycats": ["562"],
            }
        ),
        encoding="utf-8",
    )


def test_unknown_section_is_rejected():
    with pytest.raises(ValueError, match="Unknown multiqc sections"):
        effective_settings([{"sections": ["benchmark", "not_a_section"]}])


def test_cli_overrides_run_config(tmp_path):
    hydra = tmp_path / ".hydra"
    hydra.mkdir()
    (hydra / "config.yaml").write_text(
        yaml.safe_dump({"multiqc": {"sections": ["benchmark", "mag"], "stages": "all"}}),
        encoding="utf-8",
    )
    cli, extra = parse_cli(["--exclude", "mag", "--stages", "initial", "--", "--export"])
    settings = resolve_settings(tmp_path, cli)
    assert extra == ["--export"]
    assert settings["sections"] == {"benchmark"}
    assert settings["stages"] == {"initial"}


def test_rebuild_drops_disabled_section_and_keeps_upstream(tmp_path):
    _scores(tmp_path)
    mag = tmp_path / "initial_annotations_plots" / "mag_initial_ids_mqc.json"
    mag.write_text(
        json.dumps({"plot_type": "bargraph", "section_name": "Raw — MAG_ID read shares", "data": {"s": {"M1": 1}}}),
        encoding="utf-8",
    )
    table = tmp_path / "initial_annotations" / "sample.annotation.csv"
    table.parent.mkdir()
    table.write_text("seq,taxID_kraken2_0\na,562\n", encoding="utf-8")
    marker = tmp_path / ".log" / "checkpoints" / "viz_initial.done"
    marker.parent.mkdir(parents=True)
    marker.write_text("done\n", encoding="utf-8")
    before = (table.stat().st_mtime_ns, marker.stat().st_mtime_ns, table.read_text())

    full = bundle_multiqc(tmp_path)
    assert any("mag_" in path.name for path in full.glob("*_mqc.json"))
    assert (full / "report_settings.json").is_file()
    assert "MultiQC sections" in (full / "00_run_options_mqc.html").read_text()

    subset = bundle_multiqc(tmp_path, effective_settings([{"exclude": ["mag"]}]))
    assert not any("mag_" in path.name for path in subset.glob("*"))
    assert any("quality_scores" in path.name for path in subset.glob("*_mqc.json"))
    text = (subset / "98_report_status_mqc.html").read_text()
    assert "mag" not in text
    assert "composition" in text

    narrowed = bundle_multiqc(
        tmp_path,
        effective_settings([{"sections": ["benchmark"], "annotators": ["kaiju"]}]),
    )
    scores = json.loads(next(path for path in narrowed.glob("*quality_scores*") if path.suffix == ".json").read_text())
    assert "kaiju" in scores["data"]
    assert "kraken2" not in scores["data"]
    assert not any("mag_" in path.name for path in narrowed.glob("*"))
    assert (table.stat().st_mtime_ns, marker.stat().st_mtime_ns, table.read_text()) == before


def test_profile_only_mixed_and_no_ground_truth(tmp_path):
    import pandas as pd
    from samovar.profiles import write_profiles

    table = pd.DataFrame({"taxid": ["562", "0"], "N_sample": [3, 1]})
    write_profiles(
        tmp_path / "initial_profiles",
        {"baseline": table},
        {"baseline": {"normalization": "all", "database": "fixture"}},
    )
    profile_only = bundle_multiqc(tmp_path, effective_settings([{"sections": ["profiles", "composition"]}]))
    assert list(profile_only.glob("*profile_inventory*"))
    status = (profile_only / "98_report_status_mqc.html").read_text()
    assert "composition" in status
    assert "no gold standard" not in status.lower()

    _scores(tmp_path)
    mixed = bundle_multiqc(tmp_path)
    assert list(mixed.glob("*profile_inventory*"))
    assert any("quality_scores" in path.name for path in mixed.glob("*_mqc.json"))

    no_truth = bundle_multiqc(tmp_path, effective_settings([{"stages": ["initial"]}]))
    assert "regenerated" not in (no_truth / "report_settings.json").read_text().split("stages")[1].split("]")[0]
