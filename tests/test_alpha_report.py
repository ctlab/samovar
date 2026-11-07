"""Alpha diversity uses the saved profile contract and the existing Shannon."""

import json
import math

import pandas as pd

from samovar.multiqc_report import effective_settings
from samovar.profiles import alpha_metrics, metadata, write_profiles
from samovar.stage_report import bundle_multiqc
from samovar.table_scorers import shannon_entropy


def test_richness_and_shannon_exclude_unclassified():
    table = pd.DataFrame({"taxid": ["562", "9606", "0", "unmapped:native"], "N_s1": [3, 1, 4, 2]})
    info = metadata({"rank": "species", "unit": "counts", "basis": "reads", "normalization": "all"})
    rows = alpha_metrics(table, info)
    assert rows[0]["sample"] == "s1"
    assert rows[0]["richness"] == 2
    assert abs(rows[0]["shannon"] - shannon_entropy([3, 1])) < 1e-12
    assert abs(rows[0]["shannon"] - -(0.75 * math.log(0.75) + 0.25 * math.log(0.25))) < 1e-12
    assert abs(rows[0]["unclassified"] - 0.4) < 1e-12
    assert abs(rows[0]["unresolved"] - 0.2) < 1e-12


def test_zero_and_all_unclassified_and_presence():
    blank = pd.DataFrame({"taxid": ["562", "0"], "N_only": [0.0, 0.0]})
    info = metadata({"unit": "counts", "normalization": "all", "rank": "species"})
    zero = alpha_metrics(blank, info)[0]
    assert zero["richness"] == 0 and zero["shannon"] == 0.0
    unclassified = pd.DataFrame({"taxid": ["0"], "N_only": [5.0]})
    only = alpha_metrics(unclassified, info)[0]
    assert only["richness"] == 0 and only["shannon"] == 0.0
    assert abs(only["unclassified"] - 1.0) < 1e-12
    presence = metadata({"kind": "presence", "unit": "binary", "rank": "species"})
    assert alpha_metrics(pd.DataFrame({"taxid": ["562"], "N_only": [1]}), presence) is None


def test_bundle_alpha_for_several_profiles(tmp_path):
    organisms = metadata({"rank": "species", "unit": "percent", "basis": "organisms", "normalization": "all"})
    reads = metadata({"rank": "species", "unit": "counts", "basis": "reads", "normalization": "all"})
    table = pd.DataFrame({"taxid": ["562", "0"], "N_s1": [75, 25]})
    write_profiles(tmp_path / "initial_profiles", {"mpa_a": table, "mpa_b": table}, {"mpa_a": organisms, "mpa_b": reads})
    write_profiles(
        tmp_path / "reprofiled_profiles",
        {"mpa_a": table},
        {"mpa_a": {**organisms, "method": "profile_ratio"}},
    )
    write_profiles(
        tmp_path / "regenerated_profiles",
        {
            "blank": pd.DataFrame({"taxid": ["562"], "N_s1": [0.0]}),
            "flags": pd.DataFrame({"taxid": ["562"], "N_s1": [1]}),
        },
        {"blank": reads, "flags": metadata({"kind": "presence", "rank": "species"})},
    )
    staged = bundle_multiqc(tmp_path, effective_settings([{"sections": ["alpha"]}]))
    tables = list(staged.glob("*alpha*table_mqc.json"))
    names = {path.name for path in tables}
    assert any("initial_alpha_mpa_a" in name for name in names)
    assert any("initial_alpha_mpa_b" in name for name in names)
    assert any("reprofiled_alpha_mpa_a" in name for name in names)
    initial = json.loads(next(path for path in tables if "initial_alpha_mpa_a" in path.name).read_text())
    assert "natural log" in initial["description"]
    assert "basis organisms" in initial["description"]
    assert "method profile_ratio" in json.loads(
        next(path for path in tables if "reprofiled_alpha" in path.name).read_text()
    )["description"]
    row = next(iter(initial["data"].values()))
    assert row["sample"] == "s1" and row["richness"] == 1
    assert abs(row["shannon"] - shannon_entropy([0.75])) < 1e-6
    assert abs(row["unclassified"] - 0.25) < 1e-6
    assert row["basis"] == "organisms"
    other = json.loads(next(path for path in tables if "initial_alpha_mpa_b" in path.name).read_text())
    assert next(iter(other["data"].values()))["basis"] == "reads"
    csv = pd.read_csv(staged / "alpha" / "initial" / "mpa_a.csv")
    assert abs(float(csv.loc[csv["sample"].eq("s1"), "unclassified"].iloc[0]) - 0.25) < 1e-9
    blank = json.loads(next(path for path in tables if "regenerated_alpha_blank" in path.name).read_text())
    assert next(iter(blank["data"].values()))["richness"] == 0
    presence = json.loads(next(path for path in tables if "regenerated_alpha_flags" in path.name).read_text())
    assert next(iter(presence["data"].values()))["status"] == "presence"

    off = bundle_multiqc(tmp_path, effective_settings([{"sections": ["profiles"], "exclude": ["alpha"]}]))
    assert not list(off.glob("*alpha*"))
    assert list(off.glob("*profile_inventory*"))
