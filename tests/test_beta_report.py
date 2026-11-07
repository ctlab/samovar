"""Beta diversity is one full-profile distance matrix and PCoA per profile."""

import json

import numpy as np
import pandas as pd

from samovar.multiqc_report import effective_settings
from samovar.profiles import beta_diversity, classical_pcoa, grouped_beta, metadata, write_profiles
from samovar.stage_report import bundle_multiqc


def _info(**extra):
    return metadata({"rank": "species", "unit": "fraction", "basis": "reads", "normalization": "all", **extra})


def test_bray_curtis_matrix_and_pcoa_keep_negative_eigenvalues():
    table = pd.DataFrame({
        "taxid": ["1", "2", "3"],
        "N_a": [1.0, 0.0, 0.0],
        "N_b": [0.0, 1.0, 0.0],
        "N_c": [0.0, 0.0, 1.0],
        "N_d": [0.5, 0.5, 0.0],
        "N_empty": [0.0, 0.0, 0.0],
    })
    result = beta_diversity(table, _info())
    assert result["sample_ids"] == ["a", "b", "c", "d"]
    assert result["omitted_sample_ids"] == ["empty"]
    assert result["matrix"] == [
        [0.0, 1.0, 1.0, 0.5],
        [1.0, 0.0, 1.0, 0.5],
        [1.0, 1.0, 0.0, 1.0],
        [0.5, 0.5, 1.0, 0.0],
    ]
    assert len(result["explained_variance"]) >= 2
    assert abs(sum(result["explained_variance"]) - 1) < 1e-9
    assert result["status"] is None
    bent = classical_pcoa([[0, 1, 1], [1, 0, 4], [1, 4, 0]])
    assert bent["negative_eigenvalues"]
    assert all(value < 0 for value in bent["negative_eigenvalues"])
    assert 0.0 not in bent["negative_eigenvalues"]
    assert len(bent["coordinates"]) == len(bent["explained_variance"])


def test_short_empty_identical_and_presence_profiles_have_status():
    info = _info(unit="counts")
    one = beta_diversity(pd.DataFrame({"taxid": ["1"], "N_only": [1.0]}), info)
    assert one["status"] == "fewer than two samples with positive mass"
    assert one["matrix"] is None
    same = beta_diversity(pd.DataFrame({"taxid": ["1", "0"], "N_a": [1.0, 1.0], "N_b": [1.0, 1.0]}), info)
    assert same["status"] == "degenerate distances"
    assert same["matrix"] == [[0.0, 0.0], [0.0, 0.0]]
    unclassified = beta_diversity(pd.DataFrame({"taxid": ["0"], "N_a": [4.0], "N_b": [9.0]}), info)
    assert unclassified["status"] == "degenerate distances"
    assert unclassified["sample_ids"] == ["a", "b"]
    presence = beta_diversity(
        pd.DataFrame({"taxid": ["1"], "N_a": [1], "N_b": [1]}),
        metadata({"kind": "presence", "rank": "species"}),
    )
    assert presence["status"] == "presence" and presence["matrix"] is None


def test_bundle_beta_keeps_profiles_separate(tmp_path):
    organisms = _info(basis="organisms")
    reads = _info(basis="reads")
    table = pd.DataFrame({
        "taxid": ["1", "2", "3"],
        "N_a": [1.0, 0.0, 0.0],
        "N_b": [0.0, 1.0, 0.0],
        "N_c": [0.0, 0.0, 1.0],
        "N_d": [0.5, 0.5, 0.0],
    })
    write_profiles(
        tmp_path / "initial_profiles",
        {"mpa_a": table, "mpa_b": table},
        {"mpa_a": organisms, "mpa_b": reads},
    )
    write_profiles(
        tmp_path / "reprofiled_profiles",
        {"mpa_a": table},
        {"mpa_a": {**organisms, "method": "profile_ratio"}},
    )
    staged = bundle_multiqc(tmp_path, effective_settings([{"sections": ["beta"]}]))
    heat = json.loads(next(staged.glob("*initial_beta_mpa_a_heatmap_mqc.json")).read_text())
    assert heat["data"][0][1] == 1.0
    assert "basis organisms" in heat["description"]
    assert "bray_curtis" in heat["description"]
    assert "PC1" in heat["description"]
    scatter = json.loads(next(staged.glob("*initial_beta_mpa_a_pcoa_mqc.json")).read_text())
    assert set(scatter["data"]) == {"a", "b", "c", "d"}
    other = json.loads((staged / "beta" / "initial" / "mpa_b.json").read_text())
    assert other["basis"] == "reads"
    assert other["sample_ids"] == ["a", "b", "c", "d"]
    saved = json.loads((staged / "beta" / "initial" / "mpa_a.json").read_text())
    assert saved["matrix"][0][1] == 1.0
    assert "profile_ratio" in (staged / "beta" / "reprofiled" / "mpa_a.json").read_text()
    grouped = json.loads((staged / "beta" / "initial" / "grouped_0.json").read_text())
    assert grouped["profiles"] == ["mpa_a"]
    assert "mpa_b" not in grouped["profiles"]
    off = bundle_multiqc(tmp_path, effective_settings([{"sections": ["profiles"], "exclude": ["beta"]}]))
    assert not list(off.glob("*beta*"))
    assert list(off.glob("*profile_inventory*"))
    assert not np.isnan(saved["explained_variance"]).any()


def test_grouped_pcoa_is_one_space_and_names_missing_samples(tmp_path):
    organisms = _info(basis="organisms", tool="metaphlan", run_name="mpa", database="toy")
    left = pd.DataFrame({
        "taxid": ["1", "2", "3"],
        "N_a": [1.0, 0.0, 0.0],
        "N_b": [0.0, 1.0, 0.0],
        "N_c": [0.0, 0.0, 1.0],
    })
    right = pd.DataFrame({
        "taxid": ["1", "2", "3"],
        "N_a": [0.5, 0.5, 0.0],
        "N_b": [0.0, 0.0, 1.0],
    })
    genus = pd.DataFrame({"taxid": ["1", "2"], "N_a": [1.0, 0.0], "N_b": [0.0, 1.0], "N_c": [0.5, 0.5]})
    result = grouped_beta(
        {"mpa_a": left, "mpa_b": right, "genus": genus},
        {"mpa_a": organisms, "mpa_b": organisms, "genus": _info(rank="genus", basis="organisms")},
    )
    species = next(space for space in result["spaces"] if space["rank"] == "species")
    assert species["profiles"] == ["mpa_a", "mpa_b"]
    assert species["status"] is None
    labels = [(point["sample"], point["profile"]) for point in species["points"]]
    assert ("c", "mpa_a") not in labels and ("c", "mpa_b") not in labels
    assert species["omitted_samples"] == [
        {"profile": "mpa_a", "sample": "c", "reason": "not in every compatible profile"}
    ]
    left_index = labels.index(("a", "mpa_a"))
    right_index = labels.index(("a", "mpa_b"))
    assert species["matrix"][left_index][right_index] == 0.5
    assert {point["profile"] for point in species["points"]} == {"mpa_a", "mpa_b"}
    assert next(space for space in result["spaces"] if space["rank"] == "genus")["profiles"] == ["genus"]
    same = pd.DataFrame({"taxid": ["1"], "N_a": [1.0], "N_b": [1.0]})
    degenerate = grouped_beta({"mpa_a": same, "mpa_b": same}, {"mpa_a": organisms, "mpa_b": organisms})
    assert degenerate["spaces"][0]["status"] == "degenerate distances"
    write_profiles(tmp_path / "initial_profiles", {"mpa_a": left, "mpa_b": right}, {"mpa_a": organisms, "mpa_b": organisms})
    staged = bundle_multiqc(tmp_path, effective_settings([{"sections": ["beta"]}]))
    own = json.loads(next(staged.glob("*initial_beta_mpa_a_pcoa_mqc.json")).read_text())
    shared = json.loads(next(staged.glob("*initial_beta_grouped_pcoa_mqc.json")).read_text())
    assert set(own["data"]) == {"a", "b", "c"}
    assert set(shared["data"]) == {"mpa_a", "mpa_b"}
    assert set(shared["data"]["mpa_a"]) == {"a | mpa_a", "b | mpa_a"}
    summary = json.loads(next(staged.glob("*initial_beta_grouped_summary_mqc.json")).read_text())
    omitted = [row for row in summary["data"].values() if row["status"] != "ok"]
    assert omitted[0]["sample"] == "c" and omitted[0]["profile"] == "mpa_a"
    saved = json.loads((staged / "beta" / "initial" / "grouped.json").read_text())
    assert saved["points"][0]["tool"] == "metaphlan"
    assert "database" in saved["points"][0]
    off = bundle_multiqc(tmp_path, effective_settings([{"exclude": ["beta"]}]))
    assert not list(off.glob("*beta*"))
    assert not list(off.glob("beta/*"))
