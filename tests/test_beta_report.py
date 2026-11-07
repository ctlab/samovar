"""Beta diversity is one full-profile distance matrix and PCoA per profile."""

import json

import numpy as np
import pandas as pd

from samovar.multiqc_report import effective_settings
from samovar.profiles import beta_diversity, classical_pcoa, metadata, write_profiles
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
    off = bundle_multiqc(tmp_path, effective_settings([{"sections": ["profiles"], "exclude": ["beta"]}]))
    assert not list(off.glob("*beta*"))
    assert list(off.glob("*profile_inventory*"))
    assert not np.isnan(saved["explained_variance"]).any()
