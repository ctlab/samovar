"""mode=benchmark keeps the existing initial window and drops the rest."""

from pathlib import Path

import pytest


def _args(tmp_path, **extra):
    (tmp_path / "reads").mkdir(exist_ok=True)
    fields = {
        "input_config": None,
        "input_dir": str(tmp_path / "reads"),
        "output_dir": str(tmp_path / "out"),
        "threads": None,
        "cores": 1,
        "mem_mb": None,
        "custom_cli_flags": {},
    }
    fields.update(extra)
    return type("Args", (), fields)()


def test_benchmark_stops_after_initial(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    config = PipelineConfig.from_args(_args(tmp_path, mode="benchmark"))
    config.export_formats = ["abundance"]
    assert config.endpoint == "viz_initial"
    script = Path(config.generate_pipeline(str(tmp_path / "out"))).read_text()
    window = (tmp_path / "out" / ".log" / "window.env").read_text()
    assert "SAMOVAR_END=viz_initial" in window
    assert "SAMOVAR_MODE=benchmark" in window
    assert "ckpt_skip viz_initial" in script
    assert "samovar.stage_report overview" in script
    assert "ckpt_skip regenerate_reads" not in script
    assert "ckpt_skip reprofile" not in script
    assert "benchmark: export corrector logistic" in script
    assert "--export identity" in script
    assert "trained_model.joblib" not in script


def test_benchmark_keeps_an_earlier_endpoint(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    config = PipelineConfig.from_args(
        _args(tmp_path, mode="benchmark", endpoint="combine_initial")
    )
    assert config.endpoint == "combine_initial"


def test_benchmark_rejects_a_later_window(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    with pytest.raises(ValueError, match="viz_initial"):
        PipelineConfig.from_args(_args(tmp_path, mode="benchmark", endpoint="viz_reprofiled"))
    with pytest.raises(ValueError, match="viz_initial"):
        PipelineConfig.from_args(_args(tmp_path, mode="benchmark", startpoint="reprofile"))


def test_full_mode_still_reaches_reprofile(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    config = PipelineConfig.from_args(_args(tmp_path))
    assert config.mode == "full"
    assert config.endpoint == "viz_reprofiled"
    script = Path(config.generate_pipeline(str(tmp_path / "out"))).read_text()
    assert "ckpt_skip viz_reprofiled" in script
    assert "ckpt_skip regenerate_reads" in script


def test_yaml_benchmark_mode(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    spec = tmp_path / "in.yaml"
    spec.write_text("mode: benchmark\n", encoding="utf-8")
    config = PipelineConfig.from_args(_args(tmp_path, input_config=str(spec)))
    assert config.mode == "benchmark"
    assert config.endpoint == "viz_initial"
    cli = PipelineConfig.from_args(
        _args(tmp_path, input_config=str(spec), mode="full")
    )
    assert cli.mode == "full"
    assert cli.endpoint == "viz_reprofiled"


def test_benchmark_profile_only_does_not_train(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    write_config({"root": str(tmp_path), "tools": {}, "custom-flags": []}, also_repo_build=False)
    config = PipelineConfig(
        output_dir=str(tmp_path / "out"),
        mode="benchmark",
        endpoint="viz_initial",
        profile_inputs={"supplied": str(tmp_path / "profile.csv")},
        export_formats=["abundance"],
    )
    script = Path(config.generate_pipeline(str(tmp_path / "out"))).read_text()
    assert "supplied sample profiles" in script
    assert "ckpt_skip regenerate_reads" not in script
    assert "trained_model.joblib" not in script
    assert "benchmark: export corrector logistic" in script
