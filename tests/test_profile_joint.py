"""The joint correction contract uses typed profiles and regenerated truth."""
import json

import numpy as np
import pandas as pd
import pytest

from samovar.abundance_correctors import run_export
from samovar.parse_annotators import Annotation
from samovar.profiles import metadata, read_profiles, write_profiles
from samovar.scores import profile_scores


def profiles(taxa, error, prefix="s"):
    cols = [f"N_{prefix}{i}" for i in range(len(taxa))]
    truth = pd.DataFrame([taxa, 1 - np.asarray(taxa)], index=["562", "4932"], columns=cols)
    tables = {}
    for name, sign in (("a", 1), ("b", -1)):
        values = np.asarray(taxa) + sign * np.asarray(error)
        table = pd.DataFrame([.95 * values, .95 * (1 - values), [.05] * len(taxa)],
                             index=["562", "4932", "9606"], columns=cols)
        tables[name] = table.rename_axis("taxid").reset_index()
    metas = {name: metadata(dict(unit="fraction", database="db_" + name)) for name in tables}
    return Annotation.from_abundance_tables(tables, metas), truth.rename_axis("taxid").reset_index()


def fit(tmp_path):
    train, gold = profiles([.1, .3, .7, .9], [.05, -.12, .18, -.05], "train")
    target, gt = profiles([.25, .65], [.15, -.2], "heldout")
    path = tmp_path / "model.json"
    run_export("profile_joint", target, tmp_path / "out", dict(
        reference=train, reference_truth=(gold, metadata(dict(unit="fraction"))), profile_model=str(path)))
    return target, gt, path


def test_joint_learns_complementary_profiles_and_saved_apply_is_order_independent(tmp_path):
    target, gt, path = fit(tmp_path)
    tables, metas = read_profiles(tmp_path / "out")
    result = tables["profile_joint"]
    joint = profile_scores(result, gt, metas["profile_joint"], metadata(dict(unit="fraction")))
    assert joint.l1_norm.max() < 1e-10
    for name in target.abundance_tables:
        alone = profile_scores(target.abundance_tables[name], gt, target.profile_metadata[name], metadata(dict(unit="fraction")))
        assert alone.l1_norm.min() > .2
    assert result.loc[result.taxid == "9606"].iloc[:, 1:].eq(0).all().all()
    model = json.loads(path.read_text())
    assert model["weights"]["562"] == pytest.approx([1 / 1.9, 1 / 1.9])
    assert model["training_samples"] == [f"N_train{i}" for i in range(4)]
    assert set(metas["profile_joint"]["source_profiles"]) == {"a", "b"}
    shuffled = Annotation.from_abundance_tables(
        {name: table.iloc[::-1][["taxid", *reversed(table.columns[1:])]] for name, table in reversed(list(target.abundance_tables.items()))},
        dict(reversed(list(target.profile_metadata.items()))))
    # Equivalent units may change; the model schema binds abundance semantics.
    shuffled.abundance_tables["a"].iloc[:, 1:] *= 100
    shuffled.profile_metadata["a"]["unit"] = "percent"
    before = path.stat().st_mtime_ns
    run_export("profile_joint", shuffled, tmp_path / "loaded", {"profile_model": str(path)})
    assert path.stat().st_mtime_ns == before
    pd.testing.assert_frame_equal(read_profiles(tmp_path / "loaded")[0]["profile_joint"], result)


@pytest.mark.parametrize("change, message", [
    ("missing", "identities"), ("extra", "identities"), ("database", "database"),
    ("basis", "basis"), ("taxonomy", "taxonomy"), ("rank", "rank"), ("samples", "sample identities"),
])
def test_joint_rejects_incompatible_prediction_inputs(tmp_path, change, message):
    target, _, path = fit(tmp_path)
    if change == "missing":
        del target.profile_metadata["b"]
    elif change == "extra":
        target.add_profiles({"extra": target.abundance_tables["b"]}, {"extra": target.profile_metadata["b"]})
    elif change == "samples":
        target.abundance_tables["b"] = target.abundance_tables["b"].rename(columns={"N_heldout0": "N_other"})
    else:
        target.profile_metadata["b"][change] = {"database": "wrong_db", "basis": "organisms", "taxonomy": "gtdb", "rank": "genus"}[change]
    with pytest.raises(ValueError, match=message):
        run_export("profile_joint", target, tmp_path / "bad", {"profile_model": str(path)})


def test_joint_absent_taxa_unseen_taxa_and_empty_profiles(tmp_path):
    target, _, path = fit(tmp_path)
    for name in target.abundance_tables:
        table = target.abundance_tables[name]
        table.iloc[:, 1:] *= .8
        target.abundance_tables[name] = pd.concat([table, pd.DataFrame({"taxid": ["2886930"], "N_heldout0": [.2], "N_heldout1": [.2]})], ignore_index=True)
    run_export("profile_joint", target, tmp_path / "unknown", {"profile_model": str(path)})
    result = read_profiles(tmp_path / "unknown")[0]["profile_joint"]
    assert result.query("taxid == '2886930'").iloc[:, 1:].gt(0).all().all()
    for name in target.abundance_tables:
        target.abundance_tables[name] = target.abundance_tables[name].iloc[:0]
    run_export("profile_joint", target, tmp_path / "empty", {"profile_model": str(path)})
    empty = read_profiles(tmp_path / "empty")[0]["profile_joint"]
    assert empty.empty and set(empty.columns) == {"taxid", "N_heldout0", "N_heldout1"}
    # Missing rows in one profile represent absence, not a missing profiler.
    target, _, _ = fit(tmp_path)
    target.abundance_tables["b"] = target.abundance_tables["b"].query("taxid != '562'")
    run_export("profile_joint", target, tmp_path / "absent", {"profile_model": str(path)})
    assert read_profiles(tmp_path / "absent")[0]["profile_joint"].iloc[:, 1:].sum().gt(0).all()


def test_joint_model_compatibility_and_training_errors(tmp_path):
    target, gold, path = fit(tmp_path)
    saved = json.loads(path.read_text())
    for field, value in (("version", 0), ("method", "profile_ratio"), ("features", []), ("policies", {})):
        model = {**saved, field: value}
        path.write_text(json.dumps(model))
        with pytest.raises(ValueError, match="Incompatible"):
            run_export("profile_joint", target, tmp_path / "bad", {"profile_model": str(path)})
    for bad_weights in ([-1, 1], [float("nan"), 1], [0]):
        model = {**saved, "weights": {**saved["weights"], "562": bad_weights}}
        path.write_text(json.dumps(model))
        with pytest.raises(ValueError, match="coefficients"):
            run_export("profile_joint", target, tmp_path / "bad_weights", {"profile_model": str(path)})
    with pytest.raises(ValueError, match="reference_truth"):
        run_export("profile_joint", target, tmp_path / "no_gt", {"reference": target})
    with pytest.raises(ValueError, match="Presence"):
        run_export("profile_joint", target, tmp_path / "presence", dict(reference=target,
            reference_truth=(gold.assign(**{c: 1 for c in gold.columns if c != "taxid"}), metadata(dict(kind="presence")))))


def test_joint_preserves_unknown_mass_and_reports_unrecoverable_training_taxa(tmp_path):
    train, gold = profiles([.1, .3, .7, .9], [.05, -.12, .18, -.05], "train")
    target, _ = profiles([.25, .65], [.15, -.2], "heldout")
    for ann in (train, target):
        for name in ann.profile_metadata:
            ann.profile_metadata[name]["normalization"] = "all"
    for table in target.abundance_tables.values():
        table.iloc[:, 1:] *= .7
    gold.iloc[:, 1:] *= .9
    gold.loc[len(gold)] = ["9605", .1, .1, .1, .1]
    run_export("profile_joint", target, tmp_path, dict(reference=train,
        reference_truth=(gold, metadata(dict(unit="fraction", normalization="all")))))
    tables, metas = read_profiles(tmp_path)
    result = tables["profile_joint"]
    assert result.iloc[:, 1:].sum().tolist() == pytest.approx([.7, .7])
    assert list(metas["profile_joint"]["unassigned_by_sample"].values()) == pytest.approx([.3, .3])
    assert "9605" not in set(result.taxid)
    assert json.loads((tmp_path / "profile_joint.json").read_text())["unrecoverable_taxa"] == ["9605"]


def test_pipeline_joint_correction_uses_regenerated_truth_and_preserves_initial_gt(tmp_path):
    from samovar.profile_pipeline import correct, score
    initial, initial_gt = profiles([.25, .65], [.15, -.2])
    regenerated, _ = profiles([.1, .3, .7, .9], [.05, -.12, .18, -.05])
    for stage, ann in (("initial", initial), ("regenerated", regenerated)):
        write_profiles(tmp_path / f"{stage}_profiles", ann.abundance_tables, ann.profile_metadata)
    # Supply actual regenerated read-header truth; initial GT is deliberately
    # wrong and may only affect the independent benchmark.
    reads = tmp_path / "regenerated_trimmed"
    reads.mkdir()
    for sample, n in enumerate((1, 3, 7, 9)):
        lines = [f"@r{i}|taxid:{'562' if i < n else '4932'}\nACGT\n+\nIIII\n" for i in range(10)]
        (reads / f"s{sample}.fastq").write_text("".join(lines))
    truth_path = tmp_path / "initial_truth.csv"
    initial_gt.iloc[:, 1:] = 0
    initial_gt.to_csv(truth_path, index=False)
    config = dict(profile_corrector="profile_joint", profile_truth=dict(metadata(dict(unit="fraction")), path=str(truth_path)))
    correct(tmp_path, config)
    model = json.loads((tmp_path / "reprofiled_profiles/profile_joint.json").read_text())
    assert model["weights"]["562"] == pytest.approx([1 / 1.9, 1 / 1.9])
    assert score(tmp_path, "reprofiled", config)
    assert pd.read_csv(tmp_path / "reprofiled_annotations_plots/profile_quality_scores.csv").fn.eq(0).all()
