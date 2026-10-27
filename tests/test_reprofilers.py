from pathlib import Path

import pandas as pd
import pytest

from samovar.config import PipelineConfig
from samovar.paths import write_config
from samovar.reprofilers import (
    apply_saved_reprofiler,
    require_known_reprofiler,
    resolve_reprofiler,
    run_reprofiler,
)
from samovar.reprofiling import (
    IncompatiblePreprocessingError,
    predict_taxid,
    preprocess_data,
    preprocess_for_model,
    train_models,
)
from samovar.tools_import import import_tool

TOOLS = Path(__file__).resolve().parent / "tools"
LINEAR_WRAPPER = TOOLS / "linear_wrapper.py"


def _tiny_tables():
    regenerated = pd.DataFrame({
        "seq": [f"r{i}" for i in range(12)],
        "taxid_dummy": [9606, 9606, 562, 562, 9606, 562, 9606, 562, 9606, 562, 9606, 562],
        "length": [50] * 12,
        "true": [9606, 9606, 562, 562, 9606, 562, 9606, 562, 9606, 562, 9606, 562],
    })
    initial = {
        "sample.annotation": pd.DataFrame({
            "seq": ["a", "b", "c"],
            "taxid_dummy": [9606, 562, 0],
            "length": [50, 50, 50],
        })
    }
    ground = {
        "dummy": pd.DataFrame({"taxid": [9606, 562], "N_1": [10, 8]})
    }
    return regenerated, ground, initial


def test_resolve_builtin_reprofilers():
    assert resolve_reprofiler("ensemble") == ("builtin", "ensemble")
    assert resolve_reprofiler(None) == ("builtin", "random_forest")
    assert resolve_reprofiler("rf") == ("builtin", "random_forest")
    assert resolve_reprofiler("AdaBoost") == ("builtin", "adaboost")
    assert resolve_reprofiler("linear") == ("builtin", "linear")
    assert require_known_reprofiler("random_forest") == "random_forest"


def test_prepare_cores_are_training_threads(tmp_path):
    import yaml

    config = PipelineConfig(input_dir="in", output_dir=str(tmp_path), cores=6)
    paths = config.generate_configs(str(tmp_path))
    repro = yaml.safe_load(Path(paths["reprofiling"]).read_text())
    export = yaml.safe_load(Path(paths["export"]).read_text())
    assert repro["cores"] == 6
    assert repro["reprofiler"] == "random_forest"
    assert repro["gpu"] is False
    assert export["cores"] == 6
    script = Path(config.generate_pipeline(str(tmp_path))).read_text()
    assert "--cores 6" in script
    regenerated, ground, initial = _tiny_tables()
    result = run_reprofiler(
        "rf",
        regenerated=regenerated,
        ground_truth=ground,
        initial=initial,
        output_dir=tmp_path / "out",
        config={"seed": 0, "cores": repro["cores"]},
    )
    assert result.model.n_jobs == 6


def test_prepare_gpu_flag_reaches_ml(tmp_path):
    import yaml

    config = PipelineConfig(input_dir="in", output_dir=str(tmp_path), cores=2, gpu=True)
    paths = config.generate_configs(str(tmp_path))
    repro = yaml.safe_load(Path(paths["reprofiling"]).read_text())
    assert repro["gpu"] is True
    script = Path(config.generate_pipeline(str(tmp_path))).read_text()
    assert "--gpu" in script
    assert "--cores 2" in script


def test_train_models_single_method():
    regenerated, _, _ = _tiny_tables()
    processed = preprocess_data(regenerated)
    best, models, metrics, cols = train_models(
        processed, test_size=0.25, methods="random_forest"
    )
    assert list(models) == ["RandomForest"]
    assert "RandomForest" in metrics
    assert best is not None
    assert "length" in cols


@pytest.mark.parametrize("method", ["rf", "ensemble", "linear"])
def test_saved_categorical_dictionaries(method, tmp_path):
    import joblib

    regenerated = pd.DataFrame({
        "seq": [f"r{i}" for i in range(90)],
        "feat_assembly_0_MAG_ID": ["M0001", "M0002", "M0003"] * 30,
        "feat_label": ["a", "b", "c"] * 30,
        "feat_numeric": [1, 2, 3] * 30,
        "read_type": ["illumina", "hybrid", "ont"] * 30,
        "true": [562, 543, 9606] * 30,
    })
    subset = regenerated.iloc[[2, 0, 5]].drop(columns="true")
    result = run_reprofiler(
        method, regenerated=regenerated, initial={"sample": subset},
        output_dir=tmp_path / method,
        config={"feature_importance": "none", "cores": 2 if method == "linear" else 1},
    )
    model = joblib.load(tmp_path / method / "trained_model.joblib")
    assert model.samovar_preprocessing_ == result.model.samovar_preprocessing_
    assert model.samovar_preprocessing_["version"] == 1
    processed = preprocess_for_model(subset, model)
    assert processed["feat_assembly_0_MAG_ID"].tolist() == [2, 0, 2]
    assert processed["read_type"].tolist() == [2, 1, 2]
    assert processed["feat_label"].tolist() == [2, 0, 2]
    assert processed["feat_numeric"].tolist() == [3, 1, 3]
    # Numeric-looking values and missing entries never change a fitted feature's type.
    unknown = subset.copy()
    unknown["feat_assembly_0_MAG_ID"] = ["new", None, 0]
    unknown["feat_label"] = [1, 2, 3]
    unknown["read_type"] = ["new", None, ""]
    before = joblib.hash(model.samovar_preprocessing_)
    unknown = preprocess_for_model(unknown, model)
    for column in ("feat_assembly_0_MAG_ID", "feat_label", "read_type"):
        assert unknown[column].tolist() == [-1, -1, -1]
        assert preprocess_for_model(subset.drop(columns=column), model)[column].tolist() == [-1] * 3
    assert joblib.hash(model.samovar_preprocessing_) == before
    assert predict_taxid(subset, model)["taxid_SAMOVAR"].tolist() == [9606, 562, 9606]
    for other in result.models.values():
        assert other.samovar_preprocessing_ == model.samovar_preprocessing_


def test_categorical_streaming_matches_full_frame(tmp_path, monkeypatch):
    import joblib
    import samovar.reprofilers as module
    from samovar.reprofilers import stream_profile_csv

    train = pd.DataFrame({
        "feat_assembly_0_MAG_ID": ["M0001", "M0002", "M0003"] * 30,
        "read_type": ["illumina", "hybrid", "ont"] * 30,
        "feat_label": ["001", "10", "other"] * 30,
        "true": [562, 543, 9606] * 30,
    })
    # Also preserve dictionaries when the public preprocessing helper is called first.
    model, _, _, _ = train_models(preprocess_data(train))
    path = tmp_path / "model.joblib"
    joblib.dump(model, path)
    model = joblib.load(path)
    raw = pd.DataFrame({
        "seq": [f"r{i}" for i in range(20_003)],
        "feat_assembly_0_MAG_ID": ["M0001"] * 20_000 + ["M0003"] * 3,
        "read_type": ["illumina"] * 20_000 + ["ont"] * 3,
        "feat_label": ["001"] * 20_000 + ["10", "001", None],
    })
    seen = []

    def record(frame, model_path):
        seen.extend(preprocess_for_model(frame, model_path)["feat_label"].tolist())
        return predict_taxid(frame, model_path)

    monkeypatch.setattr(module, "predict_taxid", record)
    raw.to_csv(tmp_path / "sample.csv", index=False)
    stream_profile_csv(tmp_path / "sample.csv", model, tmp_path / "result.csv")
    streamed = pd.read_csv(tmp_path / "result.csv")
    full = predict_taxid(raw, str(path))
    assert streamed["taxid_SAMOVAR"].tolist() == full["taxid_SAMOVAR"].tolist()
    assert streamed["taxid_SAMOVAR"].iloc[-3:].tolist() == [9606] * 3
    assert seen == [0] * 20_000 + [1, 0, -1]


def test_legacy_and_unsupported_preprocessing(tmp_path):
    regenerated, _, initial = _tiny_tables()
    model, _, _, _ = train_models(regenerated)
    del model.samovar_preprocessing_
    assert len(predict_taxid(initial["sample.annotation"], model)) == 3
    regenerated["feat_assembly_0_MAG_ID"] = ["M0001", "M0002"] * 6
    model, _, _, _ = train_models(regenerated)
    del model.samovar_preprocessing_
    with pytest.raises(IncompatiblePreprocessingError, match="Legacy model.*Retrain"):
        predict_taxid(regenerated, model)
    # An all-numeric chunk could still come from a categorical auxiliary feature.
    regenerated = regenerated.rename(columns={"feat_assembly_0_MAG_ID": "feat_label"})
    model, _, _, _ = train_models(regenerated)
    del model.samovar_preprocessing_
    regenerated["feat_label"] = 1
    with pytest.raises(IncompatiblePreprocessingError, match="Legacy model.*Retrain"):
        predict_taxid(regenerated, model)
    model.samovar_preprocessing_ = {"version": 999, "features": {}}
    with pytest.raises(IncompatiblePreprocessingError, match="Unsupported.*version"):
        predict_taxid(regenerated, model)


def test_run_builtin_random_forest(tmp_path):
    regenerated, ground, initial = _tiny_tables()
    dest = tmp_path / "out"
    result = run_reprofiler(
        "rf",
        regenerated=regenerated,
        ground_truth=ground,
        initial=initial,
        output_dir=dest,
        config={"seed": 0},
    )
    assert (dest / "trained_model.joblib").is_file()
    assert any(dest.glob("*_reprofiled.csv"))
    table = next(iter(result.tables.values()))
    assert "taxid_SAMOVAR" in table.columns


def test_run_builtin_linear(tmp_path):
    regenerated, ground, initial = _tiny_tables()
    dest = tmp_path / "out"
    result = run_reprofiler(
        "linear",
        regenerated=regenerated,
        ground_truth=ground,
        initial=initial,
        output_dir=dest,
        config={"seed": 0, "reprofiler_flags": "--C 1 --max-iter 200"},
    )
    assert result.model is not None
    assert (dest / "trained_model.joblib").is_file()
    csvs = list(dest.glob("*_reprofiled.csv"))
    assert csvs
    assert "taxid_SAMOVAR" in pd.read_csv(csvs[0]).columns


def test_apply_saved_reprofiler_does_not_refit(tmp_path):
    regenerated, ground, initial = _tiny_tables()
    trained = tmp_path / "trained"
    run_reprofiler(
        "linear",
        regenerated=regenerated,
        ground_truth=ground,
        initial=initial,
        output_dir=trained,
        config={"seed": 0},
    )
    model = trained / "trained_model.joblib"
    mtime = model.stat().st_mtime
    initial_dir = tmp_path / "initial"
    initial_dir.mkdir()
    next(iter(initial.values())).to_csv(initial_dir / "sample.annotation.csv", index=False)
    dest = tmp_path / "applied"
    apply_saved_reprofiler(
        model_path=model,
        initial_dir=initial_dir,
        output_dir=dest,
    )
    assert (dest / "trained_model.joblib").is_file()
    assert model.stat().st_mtime == mtime
    assert any(dest.glob("*_reprofiled.csv"))
    assert "taxid_SAMOVAR" in pd.read_csv(next(dest.glob("*_reprofiled.csv"))).columns


def test_custom_linear_reprofiler(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="linwrap",
        tool_type="ml",
        exec_path=str(LINEAR_WRAPPER),
        flags="--max-iter 200",
        also_repo_build=False,
    )
    assert spec[3] == "reprofiler"
    assert spec[4] == "--max-iter 200"
    regenerated, ground, initial = _tiny_tables()
    dest = tmp_path / "profiled"
    result = run_reprofiler(
        "linwrap",
        regenerated=regenerated,
        ground_truth=ground,
        initial=initial,
        output_dir=dest,
        config={"seed": 1, "reprofiler_flags": "--C 1"},
    )
    assert result.model is not None
    assert (dest / "trained_model.joblib").is_file()
    csvs = list(dest.glob("*_reprofiled.csv"))
    assert csvs
    assert "taxid_SAMOVAR" in pd.read_csv(csvs[0]).columns


def test_prepare_reprofiler_flags(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="linwrap",
        tool_type="reprofiler",
        exec_path=str(LINEAR_WRAPPER),
        flags="--max-iter 300",
        also_repo_build=False,
    )
    (tmp_path / "reads").mkdir()
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(tmp_path / "reads"),
            "output_dir": str(tmp_path / "out"),
            "reprofiler": "linwrap",
            "tool_flags": [
                ["ml", "--C 0.5"],
                ["linwrap", "--use-priors"],
            ],
        },
    )()
    config = PipelineConfig.from_args(args)
    assert config.reprofiler == "linwrap"
    assert "--C 0.5" in (config.reprofiler_flags or "")
    assert "--use-priors" in (config.reprofiler_tool_flags or {}).get("linwrap", "")
    paths = config.generate_configs(str(tmp_path / "out"))
    text = Path(paths["reprofiling"]).read_text()
    assert "linwrap" in text
    script_text = Path(config.generate_pipeline(str(tmp_path / "out"))).read_text()
    assert "--reprofiler linwrap" in script_text
    assert "--ground-truth" in script_text


def test_reprofile_ignores_tax_column_suffix(tmp_path):
    import warnings

    from samovar.annotation_io import align_contract_columns
    from samovar.reprofilers import _align_ml_tables

    regen = pd.DataFrame({"seq": ["r"], "taxid_kaiju": [9606], "true": [9606]})
    initial = tmp_path / "initial"
    initial.mkdir()
    table = initial / "combined_annotation_table.csv"
    pd.DataFrame({"seq": ["a"], "taxID_kaiju_0": [9606], "true": [9606]}).to_csv(table, index=False)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out, _ = _align_ml_tables(regen, None, initial)
    assert not any("missing from the regenerated" in str(item.message) for item in caught)
    assert list(out.columns) == ["seq", "taxid_kaiju", "true"]

    left = pd.DataFrame({"taxID_kaiju_0": [1]})
    right = pd.DataFrame({"taxid_kaiju": [2]})
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        aligned = align_contract_columns(left, right)
    assert not caught
    assert list(aligned[0].columns) == ["taxID_kaiju_0"]
    assert list(aligned[1].columns) == ["taxid_kaiju"]

    feat_a = pd.DataFrame({"feat_kaiju_0_AA": ["x"]})
    feat_b = pd.DataFrame({"feat_kaiju_1_AA": ["y"]})
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        aligned = align_contract_columns(feat_a, feat_b)
    assert not caught
    assert list(aligned[0].columns) == ["feat_kaiju_0_AA"]

    pd.DataFrame({"seq": ["a"], "taxID_kraken2_0": [562], "true": [562]}).to_csv(table, index=False)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        out, _ = _align_ml_tables(regen, None, initial)
    assert any("taxid_kraken2" in str(item.message) for item in caught)
    assert "taxid_kraken2" in out.columns
