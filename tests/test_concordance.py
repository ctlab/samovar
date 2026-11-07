"""Concordance uses one scorer contract. Mantel is the default, not the only one."""

import json
from pathlib import Path

import pandas as pd

from samovar.concordance import score_concordance
from samovar.multiqc_report import effective_settings
from samovar.profiles import metadata, write_profiles
from samovar.stage_report import bundle_multiqc

CUSTOM = Path(__file__).parents[0] / "data" / "mean_abs_concordance.py"


def _matrix(values, ids):
    return pd.DataFrame(values, index=ids, columns=ids)


def test_mantel_is_stable_under_sample_order_and_seed():
    ids = ["a", "b", "c"]
    values = [[0.0, 0.2, 0.8], [0.2, 0.0, 0.5], [0.8, 0.5, 0.0]]
    left = _matrix(values, ids)
    right = left.loc[["c", "a", "b"], ["c", "a", "b"]]
    same = score_concordance(left, left, {"concordance": "mantel", "permutations": 99, "seed": 3})
    reordered = score_concordance(left, right, {"concordance": "mantel", "permutations": 99, "seed": 3})
    again = score_concordance(left, right, {"concordance": "mantel", "permutations": 99, "seed": 3})
    assert same["statistic"] == 1.0
    assert reordered["statistic"] == same["statistic"]
    assert again["pvalue"] == reordered["pvalue"]
    changed = _matrix([[0.0, 0.9, 0.1], [0.9, 0.0, 0.4], [0.1, 0.4, 0.0]], ids)
    assert score_concordance(left, changed, {"permutations": 99, "seed": 3})["statistic"] < 1
    short = score_concordance(left.loc[["a", "b"], ["a", "b"]], left.loc[["a", "b"], ["a", "b"]], {})
    assert short["status"] == "fewer than three shared samples"
    constant = _matrix([[0.0, 0.4, 0.4], [0.4, 0.0, 0.4], [0.4, 0.4, 0.0]], ids)
    assert score_concordance(constant, constant, {})["status"] == "constant distances"
    missing = score_concordance(left, left.drop(index="c").drop(columns="c"), {})
    assert missing["status"] == "fewer than three shared samples"
    assert missing["omitted_samples"] == ["c"]
    assert missing["statistic"] is None


def test_config_selects_custom_scorer_and_report(tmp_path):
    info = metadata({"rank": "species", "unit": "fraction", "basis": "reads", "normalization": "all", "tool": "metaphlan"})
    other = metadata({"rank": "species", "unit": "fraction", "basis": "reads", "normalization": "all", "tool": "kraken"})
    left = pd.DataFrame({"taxid": ["1", "2", "3"], "N_a": [1.0, 0.0, 0.0], "N_b": [0.0, 1.0, 0.0], "N_c": [0.5, 0.5, 0.0]})
    right = pd.DataFrame({"taxid": ["1", "2", "3"], "N_a": [1.0, 0.0, 0.0], "N_b": [0.0, 0.0, 1.0], "N_c": [0.0, 0.2, 0.8]})
    for phase in ("initial", "reprofiled"):
        write_profiles(tmp_path / f"{phase}_profiles", {"mpa": left, "kraken": right}, {"mpa": info, "kraken": other})
    configs = tmp_path / ".log" / "configs"
    configs.mkdir(parents=True)
    (configs / "config_init.yaml").write_text(f"concordance: {CUSTOM}\n", encoding="utf-8")
    staged = bundle_multiqc(tmp_path, effective_settings([{"sections": ["beta"], "stages": ["initial", "reprofiled"]}]))
    tables = list(staged.glob("*beta_concordance_table_mqc.json"))
    assert {path.name.split("_beta_")[0].split("_", 1)[-1] for path in tables} >= {"initial", "reprofiled"}
    initial = json.loads(next(path for path in tables if "initial_beta_concordance" in path.name).read_text())
    assert "not adjusted" in initial["description"]
    assert "not adjusted" in next(iter(initial["data"].values()))["pvalue_policy"]
    row = next(iter(initial["data"].values()))
    assert row["method"] == "mean_abs"
    assert row["pvalue"] is None
    assert row["statistic"] > 0
    heat = json.loads(next(staged.glob("*initial_beta_concordance_heatmap_mqc.json")).read_text())
    assert heat["xcats"] == ["mpa", "kraken"]
    saved = json.loads((staged / "beta" / "initial" / "concordance.json").read_text())
    assert saved["scorer"].endswith("mean_abs_concordance.py")
    assert saved["pairs"][0]["method"] == "mean_abs"
    own = list(staged.glob("*initial_beta_mpa_heatmap_mqc.json"))
    assert own

    mantel = bundle_multiqc(tmp_path, effective_settings([{"sections": ["beta"], "stages": ["initial"]}]))
    # The run config still selects the custom scorer, so clear it and rebuild from code default.
    (configs / "config_init.yaml").write_text("concordance: mantel\n", encoding="utf-8")
    mantel = bundle_multiqc(tmp_path, effective_settings([{"sections": ["beta"], "stages": ["initial"]}]))
    switched = json.loads(next(mantel.glob("*initial_beta_concordance_table_mqc.json")).read_text())
    assert next(iter(switched["data"].values()))["method"] == "mantel"
    assert next(iter(switched["data"].values()))["pvalue"] is not None
