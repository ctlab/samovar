from pathlib import Path

import pandas as pd
import pytest

from samovar.iss_config import parse_args, setup_iss_test
from samovar.metagenome_generators import (
    MissingMetagenomeGeneratorError,
    constant_abundance_frame,
    flags_apply_to_metagenome_generator,
    require_known_metagenome_generator,
    resolve_metagenome_generator,
)
from samovar.tools_import import import_tool


def test_resolve_metagenome_builtins():
    assert resolve_metagenome_generator("camisim") == ("builtin", "camisim")
    assert resolve_metagenome_generator("hybrid") == ("builtin", "hybrid")
    assert resolve_metagenome_generator("nanosim3") == ("builtin", "nanosim")
    assert resolve_metagenome_generator("constant_iss") == ("custom", "constant_iss")


def test_constant_abundance_frame_equalizes():
    df = pd.DataFrame({"taxid": ["562", "4932"], "N_1": [10, 90], "N_2": [0, 50]})
    out = constant_abundance_frame(df, n_reads=100)
    assert list(out["N_1"]) == [50, 50]
    assert list(out["N_2"]) == [50, 50]


def test_missing_custom_metagenome_generator(tmp_path, monkeypatch):
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    with pytest.raises(MissingMetagenomeGeneratorError, match="ghost_meta"):
        require_known_metagenome_generator("ghost_meta")


def test_reads_type_rejected_as_metagenome(tmp_path, monkeypatch):
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    script = tmp_path / "x.py"
    script.write_text("def generate(spec, metadata, config):\n    return []\n")
    import_tool(
        name="onlyreads",
        tool_type="reads",
        exec_path=str(script),
        also_repo_build=False,
    )
    with pytest.raises(MissingMetagenomeGeneratorError, match="reads_generator"):
        require_known_metagenome_generator("onlyreads")


def test_prepare_flags_metagenome_generator(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    script = tmp_path / "constant_iss.py"
    script.write_text("def generate(spec, metadata, config):\n    return []\n")
    import_tool(
        name="constant_iss",
        tool_type="meta",
        exec_path=str(script),
        flags="--model hiseq",
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
            "metagenome_generator": "constant_iss",
            "tool_flags": [["metagenome_generator", "--n-reads 200"]],
        },
    )()
    config = PipelineConfig.from_args(args)
    assert config.metagenome_generator == "constant_iss"
    assert "--n-reads 200" in (config.metagenome_generator_flags or "")
    text = Path(config.generate_configs(str(tmp_path / "out"))["annotation2iss"]).read_text()
    assert "constant_iss" in text
    assert "metagenome_generator_flags" in text


def test_generate_cli_metagenome_and_flags(tmp_path):
    args = parse_args(
        [
            "--genome_dir",
            str(tmp_path),
            "--output_dir",
            str(tmp_path / "run"),
            "--metagenome_generator",
            "camisim",
            "--flags",
            "metagenome_generator",
            "--log-mu 1",
        ]
    )
    assert args.metagenome_generator == "camisim"
    assert args.tool_flags == [["metagenome_generator", "--log-mu 1"]]


def test_setup_custom_metagenome_writes_generate_sh(tmp_path, monkeypatch):
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    script = tmp_path / "constant_iss.py"
    script.write_text("def generate(spec, metadata, config):\n    return []\n")
    import_tool(
        name="constant_iss",
        tool_type="meta",
        exec_path=str(script),
        flags="--model hiseq",
        also_repo_build=False,
    )
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "562.fna").write_text(">g\nACGT\n")
    args = parse_args(
        [
            "--genome_dir",
            str(genomes),
            "--output_dir",
            str(tmp_path / "run"),
            "--metagenome_generator",
            "constant_iss",
            "--flags",
            "constant_iss",
            "--n-reads 80",
        ]
    )
    result = setup_iss_test(args)
    yaml_text = Path(result["config"]).read_text()
    assert "constant_iss" in yaml_text
    assert "n-reads" in yaml_text or "model" in yaml_text
    sh = Path(result["pipeline"]).read_text()
    assert "run_generate_from_yaml" in sh


def test_flags_apply_meta():
    assert flags_apply_to_metagenome_generator("meta", "camisim")
    assert flags_apply_to_metagenome_generator("constant_iss", "constant_iss")
    assert not flags_apply_to_metagenome_generator("iss", "camisim")
