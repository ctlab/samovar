"""Per-taxon logistic abundance correction (Annotation → export)."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from samovar.abundance import n_sample_columns, normalize_abundance_table
from samovar.abundance_correctors import (
    apply_taxon_efficiency,
    export_logistic,
    fit_taxon_efficiency,
    require_known_export,
    resolve_export,
    run_export,
)
from samovar.config import PipelineConfig
from samovar.parse_annotators import Annotation
from samovar.paths import write_config
from samovar.tools_import import import_tool, main as import_main


def _ref_and_obs():
    """Recall(562)=0.8, recall(9606)=1.0 on the reference; observed 10/10."""
    reference = Annotation.from_long_table(
        pd.DataFrame(
            {
                "seq": [f"r{i}" for i in range(10)],
                "taxID_dummy_0": ["562", "562", "562", "562", "9606"] + ["9606"] * 5,
                "true": ["562"] * 5 + ["9606"] * 5,
                "sample": ["1"] * 10,
            }
        )
    )
    observed = Annotation.from_long_table(
        pd.DataFrame(
            {
                "seq": [f"o{i}" for i in range(20)],
                "taxID_dummy_0": ["562"] * 10 + ["9606"] * 10,
                "sample": ["1"] * 20,
            }
        )
    )
    return reference, observed


def test_resolve_builtin_export():
    assert resolve_export("logistic") == ("builtin", "logistic")
    assert resolve_export("none") == ("builtin", "identity")
    assert resolve_export("off") == ("builtin", "off")
    assert require_known_export("identity") == "identity"


def test_fit_recall_and_divide(tmp_path):
    reference, observed = _ref_and_obs()
    rates = fit_taxon_efficiency(reference, C=1e6, min_efficiency=0.01)
    assert "dummy" in rates
    dummy = rates["dummy"]
    assert dummy["562"] < dummy["9606"]
    assert dummy["562"] < 0.95
    tables = observed.to_abundance()
    corrected = apply_taxon_efficiency(tables, rates, min_efficiency=0.01)
    raw = normalize_abundance_table(tables["dummy"])
    fix = normalize_abundance_table(corrected["dummy"])
    col = n_sample_columns(raw)[0]
    raw_562 = float(raw.loc[raw["taxid"].astype(str) == "562", col].iloc[0])
    fix_562 = float(fix.loc[fix["taxid"].astype(str) == "562", col].iloc[0])
    raw_9606 = float(raw.loc[raw["taxid"].astype(str) == "9606", col].iloc[0])
    fix_9606 = float(fix.loc[fix["taxid"].astype(str) == "9606", col].iloc[0])
    assert raw_562 == 10
    assert fix_562 > raw_562
    assert abs(fix_9606 - raw_9606) < 1.5


def test_export_logistic_writes_abundance(tmp_path):
    reference, observed = _ref_and_obs()
    dest = tmp_path / "out"
    written = export_logistic(
        observed,
        dest,
        {
            "reference": reference,
            "to": "abundance",
            "C": 1e6,
            "min_efficiency": 0.01,
        },
    )
    assert Path(written).exists() or dest.is_dir()
    csvs = list(dest.glob("*.csv")) if dest.is_dir() else [Path(written)]
    assert csvs
    assert (dest / "taxon_efficiency.json").is_file()


def test_identity_when_no_true(tmp_path):
    observed = Annotation.from_long_table(
        pd.DataFrame(
            {
                "seq": ["a", "b"],
                "taxID_k_0": ["562", "9606"],
                "sample": ["1", "1"],
            }
        )
    )
    dest = tmp_path / "raw"
    run_export("none", observed, dest, {"to": "abundance"})
    table = pd.read_csv(next(dest.glob("*.csv")))
    col = n_sample_columns(table)[0]
    assert int(table[col].sum()) == 2


def test_run_export_off_is_noop(tmp_path):
    observed = Annotation.from_long_table(
        pd.DataFrame({"seq": ["a"], "taxID_k_0": ["562"], "sample": ["1"]})
    )
    assert run_export("off", observed, tmp_path / "x", {}) is None


def test_taxid_true_column_alias():
    reference = Annotation.from_long_table(
        pd.DataFrame(
            {
                "read": [f"r{i}" for i in range(8)],
                "taxID_kaiju": ["562"] * 6 + ["9606"] * 2,
                "taxID_true": ["562"] * 6 + ["9606"] * 2,
            }
        )
    )
    rates = fit_taxon_efficiency(reference, C=1e6)
    assert "kaiju" in rates
    assert rates["kaiju"]["562"] > 0.5


def test_prepare_default_logistic(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(reads),
            "output_dir": str(tmp_path / "out"),
            "kraken2": [["kraken2 /tmp/k2"]],
        },
    )()
    config = PipelineConfig.from_args(args)
    assert config.export_corrector == "logistic"
    assert config.export_formats == ["abundance"]
    script = Path(config.generate_pipeline(str(tmp_path / "out"))).read_text()
    assert "samovar.abundance_correctors" in script
    assert "--export logistic" in script


def test_custom_export_import(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    script = Path(__file__).resolve().parent / "tools" / "identity_export.py"
    spec = import_tool(
        name="echo_export",
        tool_type="export",
        exec_path=str(script),
        also_repo_build=False,
    )
    assert spec[3] == "export"
    _, observed = _ref_and_obs()
    dest = tmp_path / "custom"
    written = run_export("echo_export", observed, dest, {"to": "abundance"})
    assert Path(written).exists() or dest.exists()


def test_import_pytest_export_contract(tmp_path, monkeypatch):
    from samovar.paths import update_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    monkeypatch.setattr(
        "samovar.tools_import.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    good = Path(__file__).resolve().parent / "tools" / "identity_export.py"
    rc = import_main(
        [
            "-n",
            "echo_ex",
            "--type",
            "export",
            "--exec-path",
            str(good),
            "--pytest",
        ]
    )
    assert rc == 0
