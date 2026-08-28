from pathlib import Path

import pandas as pd

from samovar.config import PipelineConfig
from samovar.paths import write_config
from samovar.reprofilers import (
    require_known_reprofiler,
    resolve_reprofiler,
    run_reprofiler,
)
from samovar.reprofiling import preprocess_data, train_models
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
    assert resolve_reprofiler("rf") == ("builtin", "random_forest")
    assert resolve_reprofiler("AdaBoost") == ("builtin", "adaboost")
    assert resolve_reprofiler("linear") == ("builtin", "linear")
    assert require_known_reprofiler("random_forest") == "random_forest"


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
