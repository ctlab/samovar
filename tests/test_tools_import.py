import json
import stat
from pathlib import Path

from samovar.annotators_wrapper import CustomAnnotator, get_annotator_instance
from samovar.main_config import parse_tool_entry
from samovar.paths import write_config, update_config
from samovar.tools_import import import_tool, main as import_main


def test_import_tool_writes_config(tmp_path, monkeypatch):
    binary = tmp_path / "mykaiju"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="kaiju",
        tool_type="a",
        env="",
        exec_name="kaiju",
        exec_path=str(binary),
        also_repo_build=False,
    )
    assert spec[2] == str(binary.resolve())
    assert spec[3] == "annotator"
    assert spec[1] == "bash"
    loaded = parse_tool_entry(json.loads(cfg.read_text())["tools"]["kaiju"], "kaiju")
    assert loaded[2] == spec[2]
    assert loaded[3] == "annotator"


def test_import_cli_conda_prefix(tmp_path, monkeypatch):
    prefix = tmp_path / "env"
    exe = prefix / "bin" / "simulator.py"
    exe.parent.mkdir(parents=True)
    exe.write_text("#!/bin/sh\n")
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    monkeypatch.setattr(
        "samovar.tools_import.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    rc = import_main(
        [
            "-n",
            "nanosim",
            "--env",
            "conda",
            "--exec",
            "simulator.py",
            "--exec-path",
            str(prefix),
            "--type",
            "meta",
        ]
    )
    assert rc == 0
    row = json.loads(cfg.read_text())["tools"]["nanosim"]
    spec = parse_tool_entry(row, "nanosim")
    assert spec[0] == "conda"
    assert spec[3] == "metagenome_generator"
    assert Path(spec[2]).resolve() == prefix.resolve()


def test_import_flags_fifth_slot(tmp_path, monkeypatch):
    binary = tmp_path / "myboil.py"
    binary.write_text("def regenerate(annotation, metadata, config):\n    return {}\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="myboil",
        tool_type="table",
        exec_path=str(binary),
        flags="--log-mu 1 --foo bar",
        also_repo_build=False,
    )
    assert spec[3] == "table_reads_generator"
    assert spec[4] == "--log-mu 1 --foo bar"
    raw = json.loads(cfg.read_text())["tools"]["myboil"]
    assert len(raw) == 5
    assert raw[4] == "--log-mu 1 --foo bar"
    spec4 = import_tool(
        name="plain",
        tool_type="table",
        exec_path=str(binary),
        also_repo_build=False,
    )
    assert len(spec4) == 4
    raw4 = json.loads(cfg.read_text())["tools"]["plain"]
    assert len(raw4) == 4


def test_import_scoring_writes_inputs_slot(tmp_path, monkeypatch):
    script = tmp_path / "counts.py"
    script.write_text("def score(inputs, output_dir, config):\n    return None\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="counts",
        tool_type="viz",
        exec_path=str(script),
        also_repo_build=False,
    )
    assert spec[3] == "scoring"
    assert spec[5] == "*annotations"
    raw = json.loads(cfg.read_text())["tools"]["counts"]
    assert len(raw) == 6
    assert raw[4] == ""
    assert raw[5] == "*annotations"
    spec2 = import_tool(
        name="table_score",
        tool_type="scoring",
        exec_path=str(script),
        inputs="*annotations/combined_annotation_table.csv",
        flags="--tag t",
        also_repo_build=False,
    )
    assert spec2[4] == "--tag t"
    assert spec2[5] == "*annotations/combined_annotation_table.csv"


def test_import_table_scoring_group(tmp_path, monkeypatch):
    from samovar.main_config import normalize_tool_group

    assert normalize_tool_group("table-scoring") == "table_scoring"
    script = tmp_path / "bray_plugin.py"
    script.write_text(
        "def score_table(observed, generated, config):\n"
        "    return {'rank_value': 0.0, 'ok': True}\n"
    )
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="bray_plugin",
        tool_type="table-scoring",
        exec_path=str(script),
        also_repo_build=False,
    )
    assert spec[3] == "table_scoring"
    raw = json.loads(cfg.read_text())["tools"]["bray_plugin"]
    assert raw[3] == "table_scoring"
    assert len(raw) == 4


def test_imported_annotator_invokes_binary_not_custom_sh(tmp_path):
    script = tmp_path / "clf"
    script.write_text("#!/bin/sh\nexit 0\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    inst = get_annotator_instance(
        "clf",
        {
            "run_name": "clf",
            "type": "clf",
            "cmd": str(script),
            "db_path": "/db",
            "threads": 4,
        },
        {},
    )
    assert isinstance(inst, CustomAnnotator)
    cmd = inst.get_snakemake_shell_cmd("a_R1.fastq", "a_R2.fastq", ["/tmp/out.out"])
    assert "custom.sh" not in cmd
    assert str(script) in cmd
    assert "-i a_R1.fastq" in cmd
    assert "-d /db" in cmd
    assert "-p clf" not in cmd


def test_custom_annotator_appends_extra_flags(tmp_path):
    script = tmp_path / "clf"
    script.write_text("#!/bin/sh\nexit 0\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    inst = get_annotator_instance(
        "clf",
        {
            "run_name": "clf",
            "type": "clf",
            "cmd": str(script),
            "db_path": "/db",
            "threads": 2,
            "extra": "--confidence 0.1 --keep-tmp",
        },
        {},
    )
    cmd = inst.get_snakemake_shell_cmd("a_R1.fastq", "a_R2.fastq", ["/tmp/out.out"])
    assert cmd.rstrip().endswith("--confidence 0.1 --keep-tmp")


def test_prepare_merges_import_and_launch_annotator_flags(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config
    from samovar.tools_import import import_tool

    binary = tmp_path / "myclf.py"
    binary.write_text("#!/usr/bin/env python3\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="myclf",
        tool_type="annotator",
        exec_path=str(binary),
        flags="--confidence 0.1",
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
            "cmd_myclf-test": [["myclf /db --threads 8"]],
            "tool_flags": [
                ["myclf", "--keep-tmp"],
                ["annotator", "--global-ann"],
            ],
        },
    )()
    config = PipelineConfig.from_args(args)
    assert len(config.annotators) == 1
    extra = config.annotators[0].extra or ""
    assert "--confidence 0.1" in extra
    assert "--threads 8" in extra
    assert "--keep-tmp" in extra
    assert "--global-ann" in extra
    yaml_text = Path(config.generate_configs(str(tmp_path / "out"))["init_annotator"]).read_text()
    assert "--confidence 0.1" in yaml_text
    assert "--keep-tmp" in yaml_text


def test_import_pytest_blocks_bad_table_tool(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    monkeypatch.setattr(
        "samovar.tools_import.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    bad = tmp_path / "nope.py"
    bad.write_text("def regenerate(data, metadata, config):\n    return 1\n")
    rc = import_main(
        [
            "-n",
            "nope",
            "--type",
            "table",
            "--exec-path",
            str(bad),
            "--pytest",
        ]
    )
    assert rc != 0
    loaded = json.loads(cfg.read_text())
    assert "nope" not in (loaded.get("tools") or {})


def test_import_pytest_accepts_identity_table(tmp_path, monkeypatch):
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    monkeypatch.setattr(
        "samovar.tools_import.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    good = Path(__file__).resolve().parent / "tools" / "identity_table.py"
    rc = import_main(
        [
            "-n",
            "echo_tab",
            "--type",
            "table",
            "--exec-path",
            str(good),
            "--pytest",
        ]
    )
    assert rc == 0
    spec = parse_tool_entry(json.loads(cfg.read_text())["tools"]["echo_tab"], "echo_tab")
    assert spec[3] == "table_reads_generator"
