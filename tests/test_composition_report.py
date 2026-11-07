"""Composition sections show saved profile shares, not benchmark scores."""

import json

import pandas as pd

from samovar.multiqc_report import effective_settings
from samovar.profiles import (
    composition_color,
    composition_display,
    composition_fractions,
    metadata,
    write_profiles,
)
from samovar.stage_report import bundle_multiqc


def _table():
    taxa = [str(index) for index in range(1, 13)] + ["0", "unmapped:native"]
    columns = {"taxid": taxa, "N_a": [10] * 12 + [5, 1], "N_b": [1] * 12 + [20, 2]}
    return pd.DataFrame(columns)


def test_display_keeps_specials_and_matches_full_mass():
    info = metadata({"rank": "species", "normalization": "all", "unit": "counts"})
    full = composition_fractions(_table(), info)
    shown = composition_display(full, top_n=10)
    assert set(shown.loc[shown.taxid.isin(["0", "unmapped:native"]), "taxid"]) == {"0", "unmapped:native"}
    assert "Other" in set(shown.taxid)
    for column in ("N_a", "N_b"):
        assert abs(float(shown[column].sum()) - float(full[column].sum())) < 1e-9
    assert composition_color("562") == composition_color("562")
    assert composition_color("562") != composition_color("0")


def test_classified_fraction_keeps_unassigned_slice():
    table = pd.DataFrame({"taxid": ["562", "9606"], "N_only": [0.2, 0.8]})
    info = metadata({"unit": "fraction", "normalization": "classified", "rank": "species"})
    info["unassigned_by_sample"] = {"only": 0.25}
    full = composition_fractions(table, info)
    assert abs(float(full.loc[full.taxid.eq("0"), "N_only"].sum()) - 0.25) < 1e-9
    assert abs(float(full.N_only.sum()) - 1) < 1e-9


def test_bundle_composition_for_initial_and_reprofiled(tmp_path):
    table = pd.DataFrame({"taxid": ["562", "0"], "N_s1": [3, 1], "N_s2": [0, 4]})
    info = {"normalization": "all", "database": "toy", "method": "profile_ratio", "rank": "species", "tool": "metaphlan"}
    for phase in ("initial", "reprofiled"):
        write_profiles(tmp_path / f"{phase}_profiles", {"mpa": table}, {"mpa": info})
    empty = pd.DataFrame({"taxid": ["562"], "N_empty": [0.0]})
    write_profiles(
        tmp_path / "regenerated_profiles",
        {"blank": empty},
        {"blank": {"normalization": "all", "unit": "counts", "rank": "species"}},
    )
    staged = bundle_multiqc(tmp_path, effective_settings([{"sections": ["composition"]}]))
    bars = list(staged.glob("*composition*bars_mqc.json"))
    assert {path.name.split("_composition_")[0].split("_", 1)[-1] for path in bars} >= {"initial", "reprofiled"}
    initial = json.loads(next(path for path in bars if "initial_composition" in path.name).read_text())
    reprofiled = json.loads(next(path for path in bars if "reprofiled_composition" in path.name).read_text())
    assert initial["pconfig"]["stacking"] == "normal"
    assert "method profile_ratio" in initial["description"]
    assert "basis" in initial["description"] and "normalization" in initial["description"]
    assert initial["categories"]["562"]["color"] == reprofiled["categories"]["562"]["color"]
    assert initial["data"]["s1"]["562"] == 75.0
    assert initial["data"]["s1"]["0 (unclassified)"] == 25.0
    assert initial["data"]["s2"]["0 (unclassified)"] == 100.0
    assert "s2" in initial["data"]
    fractions = json.loads(next(staged.glob("*initial_composition_mpa_fractions_mqc.json")).read_text())
    by_tax = {row["taxid"]: row for row in fractions["data"].values()}
    assert by_tax["562"]["s1"] == 0.75
    csv = pd.read_csv(staged / "composition" / "initial" / "mpa.csv", dtype={"taxid": str})
    assert abs(float(csv.loc[csv.taxid.eq("562"), "N_s1"].iloc[0]) - 0.75) < 1e-9
    blank_bars = list(staged.glob("*regenerated_composition_blank_bars*"))
    assert not blank_bars
    assert list(staged.glob("*regenerated_composition_blank_fractions*"))

    again = bundle_multiqc(tmp_path, effective_settings([{"sections": ["profiles"], "exclude": ["composition"]}]))
    assert not list(again.glob("*composition*"))
    assert list(again.glob("*profile_inventory*"))
