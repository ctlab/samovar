"""Configured runtimes and shell arguments survive foreign PATH and spaces."""
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

import pytest

from samovar.config import PipelineConfig
from samovar.iss_config import ISSTestConfig
from samovar.paths import command_argv, resolve_executable

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def runtime(tmp_path, monkeypatch):
    folder = tmp_path / "runtime space's Ω"
    folder.mkdir()
    python = folder / "python"
    python.symlink_to(sys.executable)
    capture = tmp_path / "snake-args.json"
    snake = folder / "snakemake"
    snake.write_text(f'#!{sys.executable}\nimport json,sys\nfrom pathlib import Path\nPath({str(capture)!r}).write_text(json.dumps(sys.argv[1:]))\n')
    snake.chmod(0o755)
    foreign = tmp_path / "foreign"
    foreign.mkdir()
    called = tmp_path / "foreign-called"
    for name in ("python", "python3", "snakemake"):
        exe = foreign / name
        exe.write_text(f'#!/bin/sh\ntouch {shlex.quote(str(called))}\nexit 91\n')
        exe.chmod(0o755)
    config = folder / "install config.json"
    config.write_text(json.dumps({"root": str(REPO), "compilers": {"python": str(python)},
                                  "tools": {"snakemake": {"exec": {"path": str(snake)}}}}))
    monkeypatch.setenv("SAMOVAR_CONFIG", str(config))
    monkeypatch.setenv("PATH", str(foreign) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("PYTHON_PATH", str(foreign / "python3"))
    return config, python, snake, capture, called


def test_cli_config_python_beats_foreign_path_and_env(runtime):
    config, python, snake, capture, called = runtime
    result = subprocess.run([str(REPO / "bin/samovar"), "_runtime-env"], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert "PYTHON_PATH=" + str(python) in shlex.split(result.stdout)
    assert not called.exists()


@pytest.mark.parametrize("kind", ["generate", "prepare"])
def test_generated_scripts_keep_space_paths_as_one_argument(tmp_path, runtime, kind):
    config, python, snake, capture, called = runtime
    output = tmp_path / f"{kind} space's Ω"
    if kind == "generate":
        script = ISSTestConfig(genome_dir="", output_dir=str(output), host_genome="", n_samples=1).generate_pipeline(str(output))
    else:
        reads = tmp_path / "input space's Ω"
        reads.mkdir()
        (reads / "sample_R1.fastq").write_text("@r|taxid:562\nACGT\n+\nIIII\n")
        cfg = PipelineConfig(input_dir=str(reads), output_dir=str(output), cores=1, mem_mb=128,
                             endpoint="annotate_initial", run_multiqc=False)
        script = cfg.generate_pipeline(str(output))
    env = os.environ.copy()
    env["SAMOVAR_MEM_MB"] = "128"
    result = subprocess.run(["bash", script], env=env, text=True, capture_output=True)
    assert result.returncode == 0, result.stderr + result.stdout
    args = json.loads(capture.read_text())
    assert args[args.index("--directory") + 1] == str(output)
    configfile = args[args.index("--configfile") + 1]
    assert Path(configfile).is_file()
    assert "space's Ω" in configfile
    if kind == "prepare":
        assert (output / "initial/sample_R1.fastq").is_file()
    assert not called.exists()


def test_configured_executable_path_is_not_split(runtime, tmp_path, monkeypatch):
    config, python, snake, capture, called = runtime
    assert resolve_executable("snakemake", tool_key="snakemake") == str(snake)
    assert command_argv(str(snake)) == [str(snake)]
    assert command_argv(shlex.join([str(snake), "--help"])) == [str(snake), "--help"]


def test_classifier_uses_configured_binary_and_preserves_path_arguments(runtime, tmp_path):
    from samovar.annotators_wrapper import Kraken2Annotator

    config, python, snake, capture, called = runtime
    binary = python.parent / 'kraken2'
    binary.write_text(snake.read_text())
    binary.chmod(0o755)
    payload = json.loads(config.read_text())
    payload['tools']['kraken2'] = {'exec': {'path': str(binary)}}
    config.write_text(json.dumps(payload))
    folder = tmp_path / "classifier space's Ω"
    folder.mkdir()
    reads = folder / 'reads_R1.fastq'
    reads.write_text('@r\nACGT\n+\nIIII\n')
    database = folder / 'database files'
    annotator = Kraken2Annotator({'run_name': 'k2', 'db_path': str(database)}, {})
    outputs = annotator.get_expected_outputs('sample', str(folder))
    command = annotator.get_snakemake_shell_cmd(str(reads), '', outputs)
    result = subprocess.run(['bash', '-c', command], text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    args = json.loads(capture.read_text())
    assert args[args.index('--db') + 1] == str(database)
    assert args[args.index('--output') + 1] == outputs[1]
    assert args[-1] == str(reads)
    assert not called.exists()


def test_invalid_configured_python_does_not_use_foreign_path(runtime):
    config, python, snake, capture, called = runtime
    payload = json.loads(config.read_text())
    payload['compilers']['python'] = str(python.parent / 'missing-python')
    config.write_text(json.dumps(payload))
    result = subprocess.run([str(REPO / 'bin/samovar'), '_runtime-env'], text=True, capture_output=True)
    assert result.returncode != 0
    assert 'No executable Python in the SAMOVAR install config' in result.stderr
    assert not called.exists()


def test_prepared_internal_command_refreshes_configured_python(runtime):
    from samovar.annotators_wrapper import GcAnnotator

    config, python, snake, capture, called = runtime
    annotator = GcAnnotator({'cmd': '/old/install/python -m samovar.gc'}, {})
    assert shlex.split(annotator.cmd)[0] == str(python.resolve())


def test_direct_shell_wrapper_uses_configured_python(runtime, tmp_path):
    config, python, snake, capture, called = runtime
    reads = tmp_path / "reads space's Ω.fastq"
    reads.write_text('@r\nACGT\n+\nIIII\n')
    output = tmp_path / "out space's Ω.out"
    result = subprocess.run(['bash', str(REPO / 'src/annotators/custom.sh'),
                             '-i', str(reads), '-o', str(output), '-p', 'dummy'],
                            text=True, capture_output=True)
    assert result.returncode == 0, result.stderr
    assert output.read_text().strip() == 'r\t9606'
    assert not called.exists()
