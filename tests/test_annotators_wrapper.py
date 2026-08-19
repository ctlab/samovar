"""Tests for OOP annotators, CustomAnnotator, and the constant-taxID dummy."""

import os
from pathlib import Path

import pandas as pd
import pytest
import yaml

from samovar.annotators_wrapper import (
    ConstantTaxidAnnotator,
    CustomAnnotator,
    KaijuAnnotator,
    get_annotator_instance,
)
from samovar.config import setup_pipeline
from samovar.parse_annotators import Annotation, match_annotation, read_custom_raw

REPO = Path(__file__).resolve().parents[1]

import importlib.util

_CONSTANT_MOD = importlib.util.spec_from_file_location(
    "constant9606", REPO / "src" / "annotators" / "constant9606.py"
)
_constant9606 = importlib.util.module_from_spec(_CONSTANT_MOD)
_CONSTANT_MOD.loader.exec_module(_constant9606)
classify_fastq = _constant9606.classify_fastq


def _tiny_fastq(path: Path, n: int = 3, prefix: str = "read") -> None:
    records = []
    for i in range(n):
        records.append(f"@{prefix}{i}\nACGT\n+\nIIII\n")
    path.write_text("".join(records))


def test_factory_returns_custom_and_dummy_classes():
    custom = get_annotator_instance("centrifuge", {"run_name": "cf"}, {})
    assert isinstance(custom, CustomAnnotator)
    assert custom.tool_name == "centrifuge"
    dummy = get_annotator_instance("dummy", {"run_name": "dummy"}, {})
    assert isinstance(dummy, ConstantTaxidAnnotator)
    kaiju = get_annotator_instance("kaiju", {"run_name": "kaiju"}, {})
    assert isinstance(kaiju, KaijuAnnotator)


def test_match_annotation_custom_and_dummy():
    assert match_annotation("1_full_dummy.custom_dummy.out") == "dummy"
    assert match_annotation("s_run.custom_constant9606.out") == "constant9606"
    assert match_annotation("s_kaiju.kaiju.out") == "kaiju"
    assert match_annotation("notes.txt") is None


def test_constant9606_script_assigns_human_taxid(tmp_path):
    r1 = tmp_path / "a_R1.fastq"
    r2 = tmp_path / "a_R2.fastq"
    _tiny_fastq(r1)
    _tiny_fastq(r2)
    out = tmp_path / "dummy.out"
    n = classify_fastq(str(r1), str(out), taxid="9606", r2=str(r2))
    assert n == 3
    df = read_custom_raw(str(out))
    assert list(df.columns) == ["seq", "taxID"]
    assert set(df["taxID"].astype(str)) == {"9606"}
    assert set(df["seq"]) == {"read0", "read1", "read2"}


def test_constant_taxid_annotator_parse_and_shell(tmp_path):
    r1 = tmp_path / "s_R1.fastq"
    r2 = tmp_path / "s_R2.fastq"
    _tiny_fastq(r1, n=2)
    _tiny_fastq(r2, n=2)
    out = tmp_path / "s_dummy.custom_dummy.out"
    annotator = get_annotator_instance(
        "dummy",
        {"run_name": "dummy", "type": "dummy", "cmd": "dummy", "db_path": "."},
        {},
    )
    expected = annotator.get_expected_outputs("s", str(tmp_path))
    assert expected[0].endswith("s_dummy.custom_dummy.out")
    cmd = annotator.get_snakemake_shell_cmd(str(r1), str(r2), [str(out)])
    assert "constant9606.py" in cmd
    assert "--taxid 9606" in cmd
    rc = os.system(cmd)
    assert rc == 0
    parsed = annotator.parse_output(str(out))
    assert (parsed["taxID"].astype(str) == "9606").all()
    ann = Annotation({str(out): "dummy"})
    tax_cols = [c for c in ann.DataFrame.columns if c.startswith("taxID_")]
    assert tax_cols
    assert (ann.DataFrame[tax_cols[0]].astype(str) == "9606").all()


def test_prepare_wires_dummy_custom_type(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    _tiny_fastq(reads / "1_full_R1.fastq")
    _tiny_fastq(reads / "1_full_R2.fastq")
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(reads),
            "output_dir": str(tmp_path / "out"),
            "kraken2": None,
            "kaiju": None,
            "dummy": None,
            "cmd_dummy": [["dummy ."]],
        },
    )()
    result = setup_pipeline(args)
    init_cfg = yaml.safe_load(Path(result["configs"]["init_annotator"]).read_text())
    types = {run["type"] for run in init_cfg["run_config"]}
    assert "constant9606" in types
    run_names = {run["run_name"] for run in init_cfg["run_config"]}
    assert "dummy" in run_names


def test_snakemake_custom_dummy_rule(tmp_path):
    pytest.importorskip("snakemake")
    r1_dir = tmp_path / "initial"
    r1_dir.mkdir()
    _tiny_fastq(r1_dir / "sample_R1.fastq", n=4)
    _tiny_fastq(r1_dir / "sample_R2.fastq", n=4)
    out_dir = tmp_path / "reports"
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        yaml.dump(
            {
                "r1_dir": str(r1_dir),
                "r2_dir": str(r1_dir),
                "output_dir": str(out_dir),
                "run_config": [
                    {
                        "run_name": "dummy",
                        "type": "dummy",
                        "cmd": "dummy",
                        "db_path": ".",
                    }
                ],
            }
        )
    )
    import subprocess

    proc = subprocess.run(
        [
            "snakemake",
            "-s",
            str(REPO / "workflow" / "annotators" / "Snakefile"),
            "--configfile",
            str(cfg),
            "--cores",
            "1",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        pytest.fail(proc.stdout + "\n" + proc.stderr)
    outs = list(out_dir.glob("*.out"))
    assert outs
    df = pd.read_table(outs[0], header=None, names=["seq", "taxID"])
    assert (df["taxID"].astype(str) == "9606").all()
    assert len(df) == 4
