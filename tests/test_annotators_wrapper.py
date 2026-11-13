"""Annotator classes, custom tools, and the built-in constant-taxID classifier."""

import os
from pathlib import Path

import pandas as pd
import pytest
import yaml

from samovar.annotators_wrapper import (
    ConstantTaxidAnnotator,
    CustomAnnotator,
    KaijuAnnotator,
    Kraken2Annotator,
    get_annotator_instance,
)
from samovar.config import setup_pipeline
from samovar.parse_annotators import Annotation, match_annotation, read_custom_raw

REPO = Path(__file__).resolve().parents[1]

from samovar.baselines.constant_taxid import classify_fastq


def _tiny_fastq(path: Path, n: int = 3, prefix: str = "read") -> None:
    records = []
    for i in range(n):
        records.append(f"@{prefix}{i}\nACGT\n+\nIIII\n")
    path.write_text("".join(records))


def test_factory_returns_custom_and_constant_taxid_classes():
    """Factory maps centrifuge to CustomAnnotator and dummy aliases to ConstantTaxidAnnotator."""
    custom = get_annotator_instance("centrifuge", {"run_name": "cf"}, {})
    assert isinstance(custom, CustomAnnotator)
    assert custom.tool_name == "centrifuge"
    constant = get_annotator_instance("constant_taxid", {"run_name": "constant_taxid"}, {})
    assert isinstance(constant, ConstantTaxidAnnotator)
    dummy = get_annotator_instance("dummy", {"run_name": "dummy"}, {})
    assert isinstance(dummy, ConstantTaxidAnnotator)
    kaiju = get_annotator_instance("kaiju", {"run_name": "kaiju"}, {})
    assert isinstance(kaiju, KaijuAnnotator)
    assembly = get_annotator_instance("assembly", {"run_name": "assembly-test"}, {})
    from samovar.annotators_wrapper import AssemblyAnnotator

    assert isinstance(assembly, AssemblyAnnotator)


def test_match_annotation_custom_and_constant_taxid():
    assert match_annotation("1_full_dummy.custom_dummy.out") == "dummy"
    assert match_annotation("s_run.custom_constant9606.out") == "constant9606"
    assert match_annotation("s.custom_constant_taxid.out") == "constant_taxid"
    assert match_annotation("s_kaiju.kaiju.out") == "kaiju"
    assert match_annotation("1_full_run.assembly.out") == "assembly"
    assert match_annotation("s.kmer2.out") == "kmer2"
    assert match_annotation("s.gc.out") == "gc"
    assert match_annotation("s.kmer_encoder.out") == "kmer_encoder"
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


def test_constant9606_tolerates_blank_and_missing_mates(tmp_path):
    r1 = tmp_path / "blank_R1.fastq"
    r2 = tmp_path / "blank_R2.fastq"
    r1.write_text("\n")
    r2.write_text("\n")
    out = tmp_path / "blank.out"
    n = classify_fastq(str(r1), str(out), taxid="9606", r2=str(r2))
    assert n == 0
    assert out.read_text() == ""
    mixed = tmp_path / "ok_R1.fastq"
    _tiny_fastq(mixed, n=1, prefix="ok")
    out2 = tmp_path / "mixed.out"
    n2 = classify_fastq(str(mixed), str(out2), taxid="9606", r2=str(r2))
    assert n2 == 1
    missing = classify_fastq(str(mixed), str(tmp_path / "nomate.out"), r2=str(tmp_path / "nope.fastq"))
    assert missing == 1


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
    assert "constant_taxid.py" in cmd
    assert "--taxid 9606" in cmd
    rc = os.system(cmd)
    assert rc == 0
    parsed = annotator.parse_output(str(out))
    assert (parsed["taxID"].astype(str) == "9606").all()
    ann = Annotation({str(out): "dummy"})
    tax_cols = [c for c in ann.DataFrame.columns if c.startswith("taxID_")]
    assert tax_cols
    assert (ann.DataFrame[tax_cols[0]].astype(str) == "9606").all()


def test_prepare_wires_constant_taxid_from_dummy_flag(tmp_path):
    """``--dummy`` still prepares a run; type is the built-in constant_taxid annotator."""
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
    assert "constant_taxid" in types
    run_names = {run["run_name"] for run in init_cfg["run_config"]}
    assert "dummy" in run_names


def test_skip_empty_reads_and_parse_empty_kraken2(tmp_path):
    from samovar.annotators_wrapper import skip_empty_reads_cmd

    empty = tmp_path / "empty_R1.fastq"
    empty.write_text("")
    report = tmp_path / "s.kraken2.report"
    out = tmp_path / "s.kraken2.out"
    cmd = skip_empty_reads_cmd(str(empty), [str(report), str(out)], "false")
    assert "touch" in cmd
    assert "false" in cmd
    rc = os.system(cmd)
    assert rc == 0
    assert report.exists() and out.exists()
    annotator = Kraken2Annotator(
        {"run_name": "kraken2-test", "cmd": "kraken2", "db_path": "."},
        {},
    )
    parsed = annotator.parse_output(str(out))
    assert list(parsed.columns) == ["seq", "taxID"]
    assert parsed.empty


def test_kaiju_and_kraken2_commands_support_single_end_fastq(tmp_path):
    r1 = str(tmp_path / "ont_R1.fastq")
    r2 = str(tmp_path / "ont_R2.fastq")

    kaiju = KaijuAnnotator(
        {"run_name": "k", "cmd": "kaiju", "db_path": "/db/kaiju.fmi"},
        {},
    )
    kaiju_cmd = kaiju.get_snakemake_shell_cmd(
        r1, r2, [str(tmp_path / "k.out")]
    )
    assert "fastq_has_reads" in kaiju_cmd
    assert f"-j {r2}" in kaiju_cmd
    assert f"else {kaiju.cmd}" in kaiju_cmd

    kraken = Kraken2Annotator(
        {"run_name": "k2", "cmd": "kraken2", "db_path": "/db/k2"},
        {},
    )
    outputs = kraken.get_expected_outputs("ont", str(tmp_path))
    kraken_cmd = kraken.get_snakemake_shell_cmd(r1, r2, outputs)
    assert "fastq_has_reads" in kraken_cmd
    assert "--memory-mapping" in kraken_cmd
    assert "--memory-mapping --memory-mapping" not in kraken_cmd
    assert f"--paired {r1} {r2}" in kraken_cmd
    assert f"else {kraken.cmd}" in kraken_cmd
    se = kraken.get_snakemake_shell_cmd(r1, "", outputs)
    assert "--paired" not in se
    se_kaiju = kaiju.get_snakemake_shell_cmd(r1, "", [str(tmp_path / "k.out")])
    assert "-j " not in se_kaiju


def test_snakemake_constant_taxid_rule(tmp_path):
    """Snakemake custom_tool rule runs the constant-taxID annotator on toy reads."""
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
                "incremental_annotations": False,
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
    bench = next((out_dir / "runtime" / "samples").glob("*/*/*.tsv"))
    header = bench.read_text(encoding="utf-8").splitlines()[0].split("\t")
    assert "s" in header
    assert outs
    df = pd.read_table(outs[0], header=None, names=["seq", "taxID"])
    assert (df["taxID"].astype(str) == "9606").all()
    assert len(df) == 4


def test_builtin_rule_keeps_underscore_run_name(tmp_path):
    """rule kraken2 must not parse A_kraken2_x as sample A_kraken2."""
    pytest.importorskip("snakemake")
    import subprocess

    from samovar.paths import command_argv, resolve_executable, workflow_dir

    reads = tmp_path / "reads"
    reads.mkdir()
    _tiny_fastq(reads / "A.fastq", n=1, prefix="a")
    _tiny_fastq(reads / "sample_one.fastq", n=1, prefix="b")
    cfg = tmp_path / "config.yaml"
    cfg.write_text(
        yaml.dump(
            {
                "r1_dir": str(reads),
                "r2_dir": str(reads),
                "output_dir": str(tmp_path / "reports"),
                "cores": 1,
                "run_config": [
                    {
                        "run_name": "kraken2_x",
                        "type": "kraken2",
                        "cmd": "kraken2",
                        "db_path": str(tmp_path / "db"),
                        "threads": 1,
                    }
                ],
            }
        )
    )
    result = subprocess.run(
        [
            *command_argv(resolve_executable("snakemake", tool_key="snakemake")),
            "-s",
            str(Path(workflow_dir()) / "annotators" / "Snakefile"),
            "--configfile",
            str(cfg),
            "--cores",
            "1",
            "--dry-run",
        ],
        cwd=REPO,
        capture_output=True,
        text=True,
    )
    text = result.stdout + result.stderr
    assert result.returncode == 0, text
    assert "sample=A, run_name=kraken2_x" in text
    assert "sample=sample_one, run_name=kraken2_x" in text
    assert "sample=A_kraken2" not in text
@pytest.mark.parametrize("paired", [False, True])
def test_sylph_profile_command_uses_config_and_quoted_paths(tmp_path, monkeypatch, paired):
    import shlex
    from samovar.annotators_wrapper import ProfileAnnotator
    exe = tmp_path / "tools with spaces" / "sylph"
    exe.parent.mkdir()
    exe.touch()
    called = []
    def resolve(command, **kwargs):
        called.append(command)
        return command
    monkeypatch.setattr("samovar.paths.resolve_executable", resolve)
    db = tmp_path / "db with spaces"
    db.mkdir()
    (db / "genomes.syldb").write_bytes(b"index")
    (db / "genomes.syldb.taxonomy.tsv").write_text("Genome_file\ttaxid\ng\t562\n")
    (db / "sylph.manifest.json").write_text('{"compression": 50}')
    annotator = ProfileAnnotator(dict(type="sylph", cmd=str(exe), run_name="sylph_ref", db_path=str(db), threads=3), {})
    r1, r2 = str(tmp_path / "reads 1.fq"), str(tmp_path / "reads 2.fq") if paired else None
    outputs = annotator.get_expected_outputs("original_sample", str(tmp_path))
    lexer = shlex.shlex(annotator.get_snakemake_shell_cmd(r1, r2, outputs), posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    command = list(lexer)
    assert called and str(exe) in command and str(db / "genomes.syldb") in command
    command = command[command.index(str(exe)):]
    assert command[command.index("-t") + 1] == "3"
    assert command[command.index("-c") + 1] == "50"
    assert r1 in command and ("-2" in command) == paired
    assert command[command.index("--sample") + 1] == "original_sample"
    assert command[command.index("--format") + 1] == "sylph"

@pytest.mark.parametrize('paired', [False, True])
def test_bracken_profile_uses_kraken_contract_threads_and_quoted_paths(tmp_path, monkeypatch, paired):
    import shlex
    from samovar.annotators_wrapper import ProfileAnnotator
    db = tmp_path / 'database with spaces'
    db.mkdir()
    (db / 'database126mers.kmer_distrib').write_text('distribution')
    exe = tmp_path / 'programs with spaces' / 'bracken'
    exe.parent.mkdir(); exe.touch()
    monkeypatch.setattr('samovar.paths.resolve_executable', lambda cmd, **k: str(exe) if cmd == str(exe) else cmd)
    tool = ProfileAnnotator(dict(type='bracken', cmd=str(exe), run_name='kb', db_path=str(db), threads=3, extra='-t 0'), {})
    r1, r2 = str(tmp_path / 'read 1.fastq'), str(tmp_path / 'read 2.fastq') if paired else None
    output = tool.get_expected_outputs('original_sample', str(tmp_path))
    lexer = shlex.shlex(tool.get_snakemake_shell_cmd(r1, r2, output), posix=True, punctuation_chars=True)
    lexer.whitespace_split = True
    argv = list(lexer)
    assert argv[argv.index('--threads') + 1] == '3'
    assert argv[argv.index('--read-length') + 1] == '126'
    assert argv[argv.index('--threshold') + 1] == '0'
    assert argv[argv.index('--sample') + 1] == 'original_sample'
    assert r1 in argv and ('--paired' in argv) == paired
    assert 'samovar.profiles' in argv and 'bracken' in argv
    with pytest.raises(FileNotFoundError, match='distribution'):
        ProfileAnnotator({**tool.run_config, 'extra': '-r 150'}, {}).get_snakemake_shell_cmd(r1, r2, output)


@pytest.mark.parametrize('paired', [False, True])
def test_motus_profile_command_reuses_shared_profile_stage(tmp_path, monkeypatch, paired):
    import shlex
    from samovar.annotators_wrapper import ProfileAnnotator
    native = tmp_path/'motus'; native.touch()
    monkeypatch.setattr('samovar.paths.resolve_executable', lambda *a, **k: str(native))
    obj = get_annotator_instance('motus', dict(run_name='motus_named', type='motus', cmd=str(native),
         db_path=str(tmp_path/'marker DB'), threads=2, output_level='profile'), {})
    assert isinstance(obj, ProfileAnnotator)
    outputs = obj.get_expected_outputs('sample_one', str(tmp_path/'output with spaces'))
    cmd = obj.get_snakemake_shell_cmd(str(tmp_path/'r1 reads.fq'), str(tmp_path/'r2 reads.fq') if paired else '', outputs)
    args = shlex.split(cmd)
    assert 'motus' in args and 'ingest' in args and 'samovar.profiles' in args
    assert args[args.index('--sample')+1] == 'sample_one'
    assert ('--input-r2' in args) == paired
    assert args[args.index('--threads')+1] == '2'
    assert '--format' in args and args[args.index('--format')+1] == 'motus'
