"""MultiQC labels use the database token from the annotator flag."""

import json

import yaml

from samovar.config import _flag_database_label
from samovar.stage_report import (
    apply_database_labels,
    database_column_labels,
    write_bargraph_mqc,
    write_heatmap_mqc,
)


def test_flag_database_label_keeps_the_typed_token():
    assert _flag_database_label("toy") == "toy"
    assert _flag_database_label("/data/plusPF") == "plusPF"
    assert _flag_database_label(".") == ""


def test_column_index_follows_combiner_filename_order():
    pairs, singles = database_column_labels(
        [
            {"run_name": "plus", "type": "kraken2", "db_label": "plusPF"},
            {"run_name": "phage", "type": "kraken2", "db_label": "toy"},
            {"run_name": "asm", "type": "assembly", "db_path": "."},
        ]
    )
    assert pairs[("kraken2", "1")] == "toy"
    assert pairs[("kraken2", "2")] == "plusPF"
    assert ("assembly", "0") not in pairs
    assert "kraken2" not in singles
    assert (
        apply_database_labels("F1 heatmap — taxID_kraken2_1", pairs, singles)
        == "F1 heatmap — taxID_kraken2_toy"
    )
    assert (
        apply_database_labels("feat_assembly_0_MAG_ID", pairs, singles)
        == "feat_assembly_0_MAG_ID"
    )


def test_single_database_replaces_unsuffixed_tool_name():
    pairs, singles = database_column_labels(
        [{"run_name": "k2", "type": "kraken2", "db_path": "/data/plusPF"}]
    )
    assert pairs[("kraken2", "0")] == "plusPF"
    assert singles["kraken2"] == "plusPF"
    assert apply_database_labels("taxid_kraken2", pairs, singles) == "taxid_kraken2_plusPF"
    assert apply_database_labels("taxID_kraken2_0", pairs, singles) == "taxID_kraken2_plusPF"


def test_multiqc_sections_show_the_flag_database(tmp_path):
    root = tmp_path / "run"
    cfg = root / ".log" / "configs"
    cfg.mkdir(parents=True)
    (cfg / "config_init.yaml").write_text(
        yaml.safe_dump(
            {
                "run_config": [
                    {"run_name": "plus", "type": "kraken2", "db_label": "plusPF"},
                    {"run_name": "phage", "type": "kraken2", "db_label": "toy"},
                    {"run_name": "kj", "type": "kaiju", "db_path": "/db/nr"},
                ]
            }
        ),
        encoding="utf-8",
    )
    plots = root / "initial_plots"
    write_bargraph_mqc(
        {"F1": {"taxID_kraken2_1": 0.4, "taxID_kraken2_2": 0.2, "taxID_kaiju_0": 0.1}},
        plots / "importance_bars_mqc.json",
        section_name="Feature importance — taxID_kraken2_1",
    )
    bars = json.loads((plots / "importance_bars_mqc.json").read_text(encoding="utf-8"))
    assert bars["section_name"] == "Feature importance — taxID_kraken2_toy"
    assert set(bars["data"]["F1"]) == {"taxID_kraken2_toy", "taxID_kraken2_plusPF", "taxID_kaiju_nr"}
    assert bars["categories"]["taxID_kraken2_toy"]["name"] == "taxID_kraken2_toy"

    import pandas as pd

    frame = pd.DataFrame([[0.4]], index=["taxID_kaiju_0"], columns=["taxID_kraken2_2"])
    write_heatmap_mqc(
        frame,
        plots / "f1_mqc.json",
        section_name="F1 heatmap — taxID_kraken2_2",
        description="Scores for taxID_kaiju_0",
        xlab="Predicted",
        ylab="Truth",
    )
    heat = json.loads((plots / "f1_mqc.json").read_text(encoding="utf-8"))
    assert heat["section_name"] == "F1 heatmap — taxID_kraken2_plusPF"
    assert heat["description"] == "Scores for taxID_kaiju_nr"
    assert heat["xcats"] == ["taxID_kraken2_plusPF"]
    assert heat["ycats"] == ["taxID_kaiju_nr"]
    assert heat["id"]
