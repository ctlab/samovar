import json

import pandas as pd
import pytest

from samovar.resample import apply_resampling, resample_cat


def _cats():
    initial = pd.DataFrame(
        {
            "sample": ["s1", "s2", "s3"],
            "taxid": [562, 562, 9606],
            "feat": [1, 1, 3],
            "confidence": [0.9, 0.8, 0.7],
            "true": [562, 562, 9606],
            "taxid_SAMOVAR": [1, 1, 1],
        }
    )
    regenerated = pd.DataFrame(
        {
            "sample": ["s1", "s9"],
            "taxid": [562, 562],
            "feat": [1, 4],
            "confidence": [0.4, 0.5],
            "true": [562, 999],
            "taxid_SAMOVAR": [7, 8],
        }
    )
    return initial, regenerated


def test_two_groups_share_one_dataset_and_keep_whole_rows():
    initial, regenerated = _cats()
    before = regenerated.copy()
    result = resample_cat(
        initial,
        {"regen": regenerated},
        [
            {"name": "A", "samples": ["s1", "s2"], "dataset": "regen", "size": 2},
            {"name": "B", "samples": ["s3"], "dataset": "regen", "size": 2},
        ],
        columns=["taxid", "feat", "confidence"],
        seed=0,
    )
    table = result["table"]
    assert list(table["resampling_group"]) == ["A", "A", "B", "B"]
    assert table["draw_id"].is_unique
    assert set(table["source_row"]) <= {"regen:0", "regen:1"}
    assert regenerated.equals(before)
    assert "true" in table.columns and "taxid_SAMOVAR" in table.columns
    assert result["columns"] == ["taxid", "feat", "confidence"]
    assert result["schema_version"] == 1
    assert result["solver"] == "raking-v1"
    assert result["seed"] == 0
    assert result["native_profiler_resampled"] is False


def test_replacement_can_repeat_a_row():
    initial, regenerated = _cats()
    pool = regenerated.iloc[[0]].copy()
    result = resample_cat(
        initial,
        {"regen": pool},
        [{"name": "A", "samples": ["s1"], "dataset": "regen", "size": 3}],
        columns=["taxid", "feat"],
        seed=1,
    )
    assert list(result["table"]["source_row"]) == ["regen:0", "regen:0", "regen:0"]
    assert result["groups"][0]["residual_mismatch"] == 0


def test_one_sample_in_two_groups_is_an_error():
    initial, regenerated = _cats()
    with pytest.raises(ValueError, match="is in groups A and B"):
        resample_cat(
            initial,
            {"regen": regenerated},
            [
                {"name": "A", "samples": ["s1", "s2"], "dataset": "regen", "size": 1},
                {"name": "B", "samples": ["s2"], "dataset": "regen", "size": 1},
            ],
            columns=["feat"],
        )


def test_unknown_sample_is_an_error():
    initial, regenerated = _cats()
    with pytest.raises(ValueError, match="not in the initial CAT"):
        resample_cat(
            initial,
            {"regen": regenerated},
            [{"name": "A", "samples": ["missing"], "dataset": "regen", "size": 1}],
            columns=["feat"],
        )


def test_solver_is_deterministic_and_approximates_selected_joint_distribution():
    initial, regenerated = _cats()
    initial = pd.concat([initial.iloc[[0]]] * 100, ignore_index=True)
    regenerated = pd.concat(
        [regenerated.iloc[[1]]] * 99 + [regenerated.iloc[[0]]],
        ignore_index=True,
    )
    arguments = (
        initial,
        {"regen": regenerated},
        [{"name": "A", "samples": ["s1"], "dataset": "regen", "size": 100}],
    )
    first = resample_cat(*arguments, columns=["taxid", "feat"], seed=4)
    second = resample_cat(*arguments, columns=["taxid", "feat"], seed=4)
    pd.testing.assert_frame_equal(first["table"], second["table"])
    assert set(first["table"]["feat"]) == {1}
    assert first["groups"][0]["residual_mismatch"] == 0


def test_integer_taxids_match_and_confidence_stays_numeric():
    initial = pd.DataFrame(
        {
            "sample": ["s"] * 4,
            "taxid": [562.0] * 4,
            "taxid_kraken2_confidence": [0.10, 0.11, 0.10, 0.11],
        }
    )
    regenerated = pd.DataFrame(
        {
            "sample": ["g"] * 4,
            "taxid": [562, 562, 562, 999],
            "taxid_kraken2_confidence": [0.10, 0.90, 0.10, 0.90],
        }
    )
    taxids = resample_cat(
        initial,
        {"regen": regenerated},
        [{"name": "A", "samples": ["s"], "dataset": "regen", "size": 4}],
        columns=["taxid"],
        seed=1,
    )
    assert taxids["groups"][0]["support_gaps"] == []
    assert set(taxids["table"]["taxid"]) == {562}
    confidence = resample_cat(
        initial,
        {"regen": regenerated},
        [{"name": "A", "samples": ["s"], "dataset": "regen", "size": 4}],
        columns=["taxid_kraken2_confidence"],
        seed=1,
    )
    assert confidence["groups"][0]["support_gaps"] == []
    assert 0.1 in set(confidence["table"]["taxid_kraken2_confidence"])


def test_joint_gap_names_the_missing_combination():
    initial = pd.DataFrame({"sample": ["s"], "taxid": [1], "feat": ["nine"]})
    regenerated = pd.DataFrame(
        {"sample": ["g", "g"], "taxid": [1, 0], "feat": ["zero", "nine"]}
    )
    result = resample_cat(
        initial,
        {"regen": regenerated},
        [{"name": "A", "samples": ["s"], "dataset": "regen", "size": 1}],
        columns=["taxid", "feat"],
    )
    joint = next(gap for gap in result["groups"][0]["support_gaps"] if gap["column"] == "<joint>")
    assert joint["examples"] == ["1 | nine"]


def test_support_gaps_and_identity_validation():
    initial, regenerated = _cats()
    result = resample_cat(
        initial,
        {"regen": regenerated},
        [{"name": "A", "samples": ["s3"], "dataset": "regen", "size": 2}],
        columns=["taxid", "feat"],
    )
    assert {gap["column"] for gap in result["groups"][0]["support_gaps"]} >= {
        "taxid",
        "feat",
        "<joint>",
    }
    with pytest.raises(ValueError, match="cannot use"):
        resample_cat(
            initial,
            {"regen": regenerated},
            [{"name": "A", "samples": ["s1"], "dataset": "regen", "size": 1}],
            columns=["true"],
        )
    with pytest.raises(ValueError, match="explicit"):
        resample_cat(
            initial,
            {"regen": regenerated},
            [{"name": "A", "samples": ["s1"], "dataset": "regen"}],
            columns=["feat"],
        )
    with pytest.raises(ValueError, match="unique"):
        resample_cat(
            initial,
            {"regen": regenerated},
            [
                {"name": "A", "samples": ["s1"], "dataset": "regen", "size": 1},
                {"name": "A", "samples": ["s2"], "dataset": "regen", "size": 1},
            ],
            columns=["feat"],
        )


def test_apply_resampling_saves_selected_cat_diagnostics_and_abundance(tmp_path):
    initial, regenerated = _cats()
    selected, report = apply_resampling(
        regenerated,
        {
            sample: group.drop(columns=["sample"]).reset_index(drop=True)
            for sample, group in initial.groupby("sample")
        },
        {
            "seed": 9,
            "resampling": {
                "tool": "builtin",
                "columns": ["taxid", "feat"],
                "groups": [
                    {
                        "name": "all",
                        "samples": ["s1", "s2"],
                        "dataset": "regenerated",
                        "size": 4,
                    }
                ],
            },
        },
        tmp_path,
    )
    assert len(selected) == 4
    assert not set(("draw_id", "source_row")) & set(selected.columns)
    assert set(selected["resampling_group"]) == {"all"}
    assert report["seed"] == 9
    saved = pd.read_csv(tmp_path / "resampled_training_table.csv")
    assert len(saved) == 4 and saved["draw_id"].is_unique
    manifest = json.loads((tmp_path / "resampling.json").read_text())
    assert manifest["solver"] == "raking-v1"
    abundance = pd.read_csv(tmp_path / "resampled_training_abundance" / "annotator.csv")
    assert int(abundance.filter(like="N_").sum().sum()) == 4


def test_resampled_abundance_skips_samovar_and_directory_keeps_sample(tmp_path):
    from samovar.reprofilers import load_regenerated_table

    initial, regenerated = _cats()
    regenerated["taxid_SAMOVAR"] = [7, 8]
    directory = tmp_path / "regenerated_annotations"
    directory.mkdir()
    regenerated.to_csv(directory / "combined_annotation_table.csv", index=False)
    loaded = load_regenerated_table(directory)
    assert "sample" in loaded.columns
    apply_resampling(
        loaded,
        {
            sample: group.drop(columns=["sample"]).reset_index(drop=True)
            for sample, group in initial.groupby("sample")
        },
        {
            "seed": 0,
            "resampling": {
                "columns": ["taxid"],
                "groups": [
                    {
                        "name": "A",
                        "samples": ["s1"],
                        "dataset": "regenerated",
                        "size": 1,
                    }
                ],
            },
        },
        tmp_path / "out",
    )
    names = sorted(path.name for path in (tmp_path / "out" / "resampled_training_abundance").iterdir())
    assert "SAMOVAR.csv" not in names
    report = json.loads((tmp_path / "out" / "resampling.json").read_text())
    assert report["seed"] == 0


def test_apply_resampling_uses_an_imported_contract(tmp_path, monkeypatch):
    from samovar.paths import write_config
    from samovar.tools_import import import_tool

    config_path = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(config_path))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    adapter = tmp_path / "custom_resampler.py"
    adapter.write_text(
        "from samovar.resample import resample_cat\n"
        "def resample(initial, regenerated, groups, config):\n"
        "    result = resample_cat(initial, regenerated, groups, "
        "columns=config['columns'], seed=config['seed'])\n"
        "    result['solver'] = 'custom-test'\n"
        "    return result\n"
    )
    import_tool(
        name="custom_resampler",
        tool_type="resampling",
        exec_path=str(adapter),
        also_repo_build=False,
    )
    initial, regenerated = _cats()
    _, report = apply_resampling(
        regenerated,
        {
            sample: group.drop(columns=["sample"]).reset_index(drop=True)
            for sample, group in initial.groupby("sample")
        },
        {
            "resampling": {
                "tool": "custom_resampler",
                "columns": ["taxid"],
                "groups": [
                    {
                        "name": "one",
                        "samples": ["s1"],
                        "dataset": "regenerated",
                        "size": 1,
                    }
                ],
            }
        },
        tmp_path / "out",
    )
    assert report["solver"] == "custom-test"
