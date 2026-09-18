import json
import stat
from pathlib import Path

from samovar.annotators_wrapper import CustomAnnotator, get_annotator_instance
from samovar.main_config import parse_tool_entry
from samovar.paths import write_config, update_config
from samovar.tools_import import import_tool, main as import_main
from samovar.tool_spec import bare_tool_name


def _tool_row(tools: dict, name: str):
    if name in tools:
        return tools[name]
    for key, value in tools.items():
        if bare_tool_name(key) == name:
            return value
    raise KeyError(name)


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
    tools = json.loads(cfg.read_text())["tools"]
    loaded = parse_tool_entry(_tool_row(tools, "kaiju"), "kaiju")
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
    row = _tool_row(json.loads(cfg.read_text())["tools"], "nanosim")
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
    raw = _tool_row(json.loads(cfg.read_text())["tools"], "myboil")
    assert isinstance(raw, dict)
    assert raw["flags"] == "--log-mu 1 --foo bar"
    assert "exec" in raw
    spec4 = import_tool(
        name="plain",
        tool_type="table",
        exec_path=str(binary),
        also_repo_build=False,
    )
    assert len(spec4) == 4
    raw4 = _tool_row(json.loads(cfg.read_text())["tools"], "plain")
    assert raw4["flags"] == ""


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
    raw = _tool_row(json.loads(cfg.read_text())["tools"], "counts")
    assert raw["inputs"] == "*annotations"
    assert raw.get("flags") == ""
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


def test_import_assembly_groups(tmp_path, monkeypatch):
    from samovar.main_config import normalize_tool_group

    assert normalize_tool_group("assembler") == "assembler"
    assert normalize_tool_group("prodigal") == "gene_caller"
    assert normalize_tool_group("gene-caller") == "gene_caller"
    assert normalize_tool_group("binner-qc") == "binner_qc"
    assert normalize_tool_group("binner-combine") == "binner_combine"
    assert normalize_tool_group("mag-taxonomy") == "mag_taxonomy"
    assert normalize_tool_group("aligner") == "aligner"
    assert normalize_tool_group("mag-quantifier") == "mag_quantifier"
    assert normalize_tool_group("taxon-quantifier") == "taxon_quantifier"
    assert normalize_tool_group("read-assigner") == "read_assigner"
    script = tmp_path / "dummy_assembler.py"
    script.write_text("def main(argv=None):\n    return 0\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="megahit",
        tool_type="assembler",
        exec_path=str(script),
        also_repo_build=False,
    )
    assert spec[3] == "assembler"


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
    raw = _tool_row(json.loads(cfg.read_text())["tools"], "bray_plugin")
    assert raw["type"] == "table_scoring"


def test_import_sample_scoring_group(tmp_path, monkeypatch):
    from samovar.main_config import normalize_tool_group

    assert normalize_tool_group("sample-qc") == "sample_scoring"
    script = tmp_path / "qc_plugin.py"
    script.write_text(
        "def score_samples(generated, reference, config=None):\n"
        "    return {'s1': 1.0}\n"
    )
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="qc_plugin",
        tool_type="sample-score",
        exec_path=str(script),
        also_repo_build=False,
    )
    assert spec[3] == "sample_scoring"
    raw = _tool_row(json.loads(cfg.read_text())["tools"], "qc_plugin")
    assert raw["type"] == "sample_scoring"


def test_import_sample_filtering_group(tmp_path, monkeypatch):
    from samovar.main_config import normalize_tool_group

    assert normalize_tool_group("sample-filter") == "sample_filtering"
    script = tmp_path / "filt_plugin.py"
    script.write_text(
        "def filter_samples(table, scores, config=None):\n"
        "    return table\n"
    )
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="filt_plugin",
        tool_type="sample-filter",
        exec_path=str(script),
        also_repo_build=False,
    )
    assert spec[3] == "sample_filtering"
    raw = _tool_row(json.loads(cfg.read_text())["tools"], "filt_plugin")
    assert raw["type"] == "sample_filtering"


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


def test_import_flags_translate_and_lazy_install(tmp_path, monkeypatch):
    binary = tmp_path / "myboil.py"
    binary.write_text("def regenerate(annotation, metadata, config):\n    return {}\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="myboil",
        tool_type="table",
        exec_path=str(binary),
        flags="--log-mu 1",
        flags_translate="--threads:--n-jobs --cores:--n-jobs",
        lazy_install="pip install myboil",
        version="1.2.3",
        also_repo_build=False,
    )
    tools = json.loads(cfg.read_text())["tools"]
    assert "myboil:1.2.3" in tools
    rec = tools["myboil:1.2.3"]
    assert rec["lazy-install"] == "pip install myboil"
    assert rec["flags-translate"]["--threads"] == "--n-jobs"
    assert rec["type"] == "table_reads_generator"


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
    good = Path(__file__).resolve().parents[1] / "src" / "samovar" / "baselines" / "identity_table.py"
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
    spec = parse_tool_entry(_tool_row(json.loads(cfg.read_text())["tools"], "echo_tab"), "echo_tab")
    assert spec[3] == "table_reads_generator"


def test_import_pytest_maps_every_importable_group_to_one_contract_test():
    """``samovar tools import --pytest`` must hit exactly one contract test per group."""
    from samovar.main_config import TOOL_GROUPS
    from samovar.tool_contracts import GROUP_TO_TESTNODE, _contract_repo_root, default_tool_path

    skip = {"runtime", "compiler", "workflow"}
    missing = [g for g in TOOL_GROUPS if g not in skip and g not in GROUP_TO_TESTNODE]
    extra = [g for g in GROUP_TO_TESTNODE if g not in TOOL_GROUPS]
    assert missing == [], f"import --pytest has no contract test for {missing}"
    assert extra == [], f"GROUP_TO_TESTNODE has unknown groups {extra}"

    root = _contract_repo_root()
    src = (root / "tests" / "test_tool_contracts.py").read_text(encoding="utf-8")
    for group, node in GROUP_TO_TESTNODE.items():
        rel, func = node.split("::", 1)
        assert (root / rel).is_file(), f"{group}: missing {rel}"
        assert f"def {func}(" in src, f"{group}: {func} not defined in {rel}"
        assert default_tool_path(group).is_file(), f"{group}: baseline/example dest missing"


def test_import_pytest_runs_baseline_for_every_contract_group():
    """Each group’s baseline dest must pass the same pytest node import --pytest uses."""
    from samovar.tool_contracts import GROUP_TO_TESTNODE, default_tool_path, run_contract_pytest

    failed = []
    for group in GROUP_TO_TESTNODE:
        dest = default_tool_path(group)
        code, output = run_contract_pytest(str(dest), group)
        if code != 0:
            failed.append((group, code, output[-2500:]))
    assert not failed, failed


BIB_ONE = """@article{toyA2024,
  title={Toy annotator A},
  author={Doe, Jane},
  year={2024},
  journal={None}
}
"""

BIB_TWO = """@article{toyB2024,
  title={Toy annotator B},
  year={2024}
}
@inproceedings{toyBtalk,
  title={Talk},
  year={2023}
}
"""


def test_import_inline_bibtex_writes_registry_and_record(tmp_path, monkeypatch):
    binary = tmp_path / "myclf"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    spec = import_tool(
        name="myclf",
        tool_type="annotator",
        exec_path=str(binary),
        bibtex=[BIB_ONE],
        also_repo_build=False,
    )
    assert spec[3] == "annotator"
    rec = _tool_row(json.loads(cfg.read_text())["tools"], "myclf")
    assert rec["citation"] == ["myclf_toyA2024.bib"]
    registry = json.loads((tmp_path / "cite" / "citations.json").read_text())
    assert registry["myclf"] == ["myclf_toyA2024.bib"]
    assert "@article{toyA2024" in (tmp_path / "cite" / "myclf_toyA2024.bib").read_text()


def test_import_bibtex_file_and_idempotent_merge(tmp_path, monkeypatch):
    binary = tmp_path / "other"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    bib = tmp_path / "paper.bib"
    bib.write_text(BIB_TWO)
    extra = tmp_path / "note.txt"
    extra.write_text("@misc{note, title={Note}, year={2020}}\n")
    import_tool(
        name="other",
        tool_type="annotator",
        exec_path=str(binary),
        bibtex_file=[str(bib), str(extra)],
        also_repo_build=False,
    )
    import_tool(
        name="other",
        tool_type="annotator",
        exec_path=str(binary),
        bibtex_file=[str(bib)],
        also_repo_build=False,
    )
    rec = _tool_row(json.loads(cfg.read_text())["tools"], "other")
    assert rec["citation"] == ["other_paper.bib", "other_note.txt"]
    registry = json.loads((tmp_path / "cite" / "citations.json").read_text())
    assert registry["other"] == ["other_paper.bib", "other_note.txt"]
    linked = tmp_path / "cite" / "other_paper.bib"
    assert linked.is_file() or linked.is_symlink()
    assert "toyB2024" in linked.read_text()


def test_import_citations_do_not_drop_other_tools(tmp_path, monkeypatch):
    a = tmp_path / "a.py"
    a.write_text("#!/bin/sh\n")
    b = tmp_path / "b.py"
    b.write_text("#!/bin/sh\n")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    import_tool(
        name="alpha",
        tool_type="annotator",
        exec_path=str(a),
        bibtex=[BIB_ONE],
        also_repo_build=False,
    )
    import_tool(
        name="beta",
        tool_type="annotator",
        exec_path=str(b),
        bibtex=["@article{beta, title={B}, year={2021}}\n"],
        also_repo_build=False,
    )
    registry = json.loads((tmp_path / "cite" / "citations.json").read_text())
    assert "alpha" in registry and "beta" in registry
    assert registry["alpha"] == ["alpha_toyA2024.bib"]


def test_import_cli_bibtex_and_file(tmp_path, monkeypatch):
    binary = tmp_path / "cli_tool"
    binary.write_text("#!/bin/sh\nexit 0\n")
    binary.chmod(binary.stat().st_mode | stat.S_IEXEC)
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    monkeypatch.setattr(
        "samovar.tools_import.update_config",
        lambda updates, also_repo_build=True: update_config(updates, also_repo_build=False),
    )
    bib = tmp_path / "fromfile.bib"
    bib.write_text("@article{fromfile, title={F}, year={2019}}\n")
    rc = import_main(
        [
            "-n",
            "cli_tool",
            "--type",
            "annotator",
            "--exec-path",
            str(binary),
            "--bibtex",
            BIB_ONE,
            "--bibtex-file",
            str(bib),
        ]
    )
    assert rc == 0
    rec = _tool_row(json.loads(cfg.read_text())["tools"], "cli_tool")
    assert "cli_tool_toyA2024.bib" in rec["citation"]
    assert "cli_tool_fromfile.bib" in rec["citation"]


def test_import_without_bibtex_leaves_existing_cli(tmp_path, monkeypatch):
    """Database import still works when citation flags are unused."""
    from samovar.tools_import import import_database

    db = tmp_path / "idx"
    db.write_text("x")
    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}, "databases": {}}, also_repo_build=False)
    rec = import_database(
        name="toy",
        tool="kraken2",
        exec_path=str(db),
        also_repo_build=False,
    )
    assert rec["path"]
    assert (tmp_path / "cite" / "citations.json").is_file() is False

