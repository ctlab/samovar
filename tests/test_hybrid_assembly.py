"""Hybrid libraries survive the pipeline and reach native assembler CLIs."""
import json
import subprocess
import sys
from pathlib import Path

import pytest

from samovar import assembly_profiling as ap
from samovar.camisim import _merge_hybrid_reads
from samovar.qc import trim_directory
from samovar.seqio import (
    assembly_read_libraries, drop_r2_if_source_single_end, iter_fastq_records,
    link_or_copy_reads, list_fastq_samples,
)


def _hybrid(directory):
    directory.mkdir(parents=True, exist_ok=True)
    for tech, mate, seq in (("illumina", 1, "ACGT"), ("illumina", 2, "TGCA"), ("ont", 1, "ACGTACGT")):
        (directory / f"1_{tech}_R{mate}.fastq").write_text(
            f"@{tech}/{mate}|taxid:562|read_type:{tech}\n{seq}\n+\n{'I' * len(seq)}\n"
        )
    _merge_hybrid_reads({"illumina": [], "ont": []}, directory, 1)
    return directory / "1_full_R1.fastq"


def test_hybrid_generation_staging_qc_and_assembly(tmp_path, monkeypatch):
    source = tmp_path / "source"
    _hybrid(source)
    initial = tmp_path / "run" / "initial"
    assert link_or_copy_reads(source, initial) == 4
    trimmed = initial.with_name("initial_trimmed")
    trim_directory(initial, trimmed)
    assert list_fastq_samples(trimmed) == ["1_full"]
    assert (trimmed / "1_illumina_R2.fastq").is_file()
    calls = []

    def run(cmd):
        calls.append(cmd)
        short = Path(cmd[cmd.index("-1") + 1])
        mate = Path(cmd[cmd.index("-2") + 1])
        ont = Path(cmd[cmd.index("--nanopore") + 1])
        assert short.name == "1_illumina_R1.fastq"
        assert mate.name == "1_illumina_R2.fastq"
        assert len(list(iter_fastq_records(short))) == len(list(iter_fastq_records(mate))) == 1
        assert all("read_type:ont" in r[0] for r in iter_fastq_records(ont))
        work = Path(cmd[cmd.index("-o") + 1])
        work.mkdir(parents=True)
        (work / "contigs.fasta").write_text(">c\nACGT\n")

    monkeypatch.setattr(ap, "which_tool", lambda name: "/configured/spades.py")
    monkeypatch.setattr(ap, "_run", run)
    dest = tmp_path / "contigs.fa"
    ap.run_assembler(str(trimmed / "1_full_R1.fastq"), dest=str(dest), name="spades", extra=["--isolate"])
    assert calls[0][0] == "/configured/spades.py"
    assert "--isolate" in calls[0]
    assert dest.read_text() == ">c\nACGT\n"


@pytest.mark.parametrize("assembler,expected", [("megahit", "illumina"), ("flye", "ont"), ("spades", "both")])
def test_legacy_flattened_hybrid_routes_types(tmp_path, monkeypatch, assembler, expected):
    r1 = _hybrid(tmp_path / "reads")
    for path in r1.parent.glob("1_illumina*"):
        path.unlink()
    (r1.parent / "1_ont_R1.fastq").unlink()
    observed = []

    def run(cmd):
        read_args = ("-r",) if assembler == "megahit" else (("-s", "--nanopore") if assembler == "spades" else ("--nano-hq",))
        for flag in read_args:
            reads = Path(cmd[cmd.index(flag) + 1])
            observed.extend(r[0].split("read_type:")[1].strip() for r in iter_fastq_records(reads))
        work = Path(cmd[cmd.index("-o") + 1])
        work.mkdir(parents=True, exist_ok=True)
        name = {"spades": "contigs.fasta", "megahit": "final.contigs.fa", "flye": "assembly.fasta"}[assembler]
        (work / name).write_text(">c\nACGT\n")

    monkeypatch.setattr(ap, "which_tool", lambda name, default="": "/configured/" + name)
    monkeypatch.setattr(ap, "_run", run)
    ap.run_assembler(str(r1), dest=str(tmp_path / "contigs.fa"), name=assembler)
    assert set(observed) == ({"illumina", "ont"} if expected == "both" else {expected})
    assert len(observed) == (3 if expected == "both" else 2 if expected == "illumina" else 1)


@pytest.mark.parametrize("tech,assembler", [("ont", "megahit"), ("illumina", "flye"), ("ont", "spades")])
def test_incompatible_reads_rejected_before_launch(tmp_path, monkeypatch, tech, assembler):
    r1 = tmp_path / "sample.fastq"
    r1.write_text(f"@a|read_type:{tech}\nACGT\n+\nIIII\n")
    monkeypatch.setattr(ap, "_run", lambda cmd: pytest.fail("must reject before launching"))
    with pytest.raises(ValueError, match="requires"):
        ap.run_assembler(str(r1), dest=str(tmp_path / "out.fa"), name=assembler)


@pytest.mark.parametrize("gzip_reads", [False, True])
def test_regenerated_hybrid_retains_paired_library(tmp_path, monkeypatch, gzip_reads):
    from samovar import camisim
    from samovar.reads_generators import _simulate_camisim_sample

    def generate(config):
        _hybrid(Path(config).parent / "initial")
    monkeypatch.setattr(camisim, "run_from_config", generate)
    run = tmp_path / "run"
    _hybrid(run / "initial")
    out = run / "regenerated"
    out.mkdir()
    _simulate_camisim_sample("hybrid", {"562": 3}, {"562": str(tmp_path / "562.fa")}, str(out),
                             "1_full", "kraken2", model="hiseq", cpus=1, seed=42,
                             gzip_reads=gzip_reads, flags=[], ids=[], config={})
    drop_r2_if_source_single_end(out)
    assert list_fastq_samples(out) == ["1_full_kraken2"]
    ext = ".fastq.gz" if gzip_reads else ".fastq"
    assert (out / ("1_full_kraken2_illumina_R2" + ext)).is_file()
    with assembly_read_libraries(out / ("1_full_kraken2_R1" + ext)) as (r1, r2, ont, typed):
        assert typed and r2 and ont
        assert len(list(iter_fastq_records(r1))) == len(list(iter_fastq_records(r2))) == 1


def test_hybrid_mapping_uses_separate_presets(tmp_path, monkeypatch):
    r1 = _hybrid(tmp_path / "reads")
    calls = []
    commands = []

    def align(r1, r2, reference, dest, threads, preset, extra):
        calls.append((Path(r1).name, Path(r2).name if r2 else "", preset))
        Path(dest).write_text("bam")
    monkeypatch.setattr(ap, "_run_minimap2", align)
    monkeypatch.setattr(ap, "which_tool", lambda name: "/configured/" + name)
    monkeypatch.setattr(ap, "_run", lambda cmd: commands.append(cmd))
    ap.run_minimap2(str(r1), "", "contigs.fa", str(tmp_path / "map.bam"))
    assert calls == [("1_illumina_R1.fastq", "1_illumina_R2.fastq", "sr"), ("1_ont_R1.fastq", "", "map-ont")]
    assert commands[0][1] == "merge"
    assert commands[1][1] == "index"


def test_assembly_hybrid_shell_uses_configured_tools(tmp_path, monkeypatch):
    """Exercise the real shell/Python boundary with tiny native-CLI stand-ins."""
    r1 = _hybrid(tmp_path / "reads")
    log = tmp_path / "calls.jsonl"
    monkeypatch.setenv("HYBRID_TEST_LOG", str(log))
    script = r'''import json, os, sys
from pathlib import Path
args = sys.argv[1:]
name = Path(sys.argv[0]).name
with open(os.environ["HYBRID_TEST_LOG"], "a") as f:
    f.write(json.dumps([name] + args) + "\n")
def value(flag): return args[args.index(flag) + 1]
if name == "spades.py":
    out = Path(value("-o")); out.mkdir(parents=True)
    assert "--nanopore" in args and "-1" in args and "-2" in args
    (out / "contigs.fasta").write_text(">c\nACGT\n")
elif name == "minimap2":
    reads = Path(args[-1])
    print("@SQ\tSN:c\tLN:100")
    lines = reads.read_text().splitlines()
    for h, seq in zip(lines[::4], lines[1::4]):
        print(f"{h[1:]}\t0\tc\t1\t60\t{len(seq)}M\t*\t0\t0\t{seq}\t*")
elif name == "samtools":
    if args[0] == "sort": Path(value("-o")).write_text(sys.stdin.read())
    elif args[0] == "merge": Path(args[-3]).write_text(Path(args[-2]).read_text() + Path(args[-1]).read_text())
    elif args[0] == "view": Path(value("-o")).write_text(Path(args[-1]).read_text())
    elif args[0] == "index": Path(args[1] + ".bai").touch()
elif name == "kraken2":
    Path(value("--output")).write_text("C\tc\tBacteria (taxid 562)\t100\t562:1\n")
'''
    tools = {}
    for name in ("spades.py", "minimap2", "samtools", "kraken2"):
        exe = tmp_path / name
        exe.write_text(f"#!{sys.executable}\n" + script)
        exe.chmod(0o755)
        tools["spades" if name == "spades.py" else name] = {"exec": {"path": str(exe)}}
    config = tmp_path / "config.json"
    config.write_text(json.dumps({"tools": tools}))
    monkeypatch.setenv("SAMOVAR_CONFIG", str(config))
    monkeypatch.setenv("PYTHON_PATH", sys.executable)
    out = tmp_path / "hybrid.out"
    subprocess.run(["bash", str(Path(ap.__file__).parents[1] / "annotators" / "custom.sh"),
                    "-i", str(r1), "-d", str(tmp_path / "db"), "-o", str(out),
                    "-p", "assembly_hybrid", "-t", "2"], check=True, capture_output=True, text=True)
    assert "\t562\n" in out.read_text()
    calls = [json.loads(line) for line in log.read_text().splitlines()]
    assert calls[0][0] == "spades.py"
    assert [c[c.index("-ax") + 1] for c in calls if c[0] == "minimap2"] == ["sr", "map-ont"]


@pytest.mark.optional
def test_hybrid_mapping_with_installed_minimap2_and_samtools(tmp_path):
    try:
        ap.which_tool("minimap2")
        samtools = ap.which_tool("samtools")
    except FileNotFoundError as exc:
        pytest.skip(str(exc))
    r1 = _hybrid(tmp_path / "reads")
    reference = tmp_path / "contigs.fa"
    reference.write_text(">c\n" + "ACGT" * 100 + "\n")
    bam = tmp_path / "map.bam"
    ap.run_minimap2(str(r1), "", str(reference), str(bam), threads=1)
    records = subprocess.check_output([samtools, "view", str(bam)], text=True).splitlines()
    assert len(records) == 3
    assert sum("read_type:illumina" in rec for rec in records) == 2
    assert sum("read_type:ont" in rec for rec in records) == 1
    assert bam.with_suffix(".bam.bai").is_file()


def test_partial_companions_do_not_drop_reads(tmp_path):
    r1 = _hybrid(tmp_path / "reads")
    (r1.parent / "1_ont_R1.fastq").unlink()
    with assembly_read_libraries(r1) as (short, mate, ont, typed):
        assert typed and not mate
        assert len(list(iter_fastq_records(short))) == 2
        assert len(list(iter_fastq_records(ont))) == 1


def test_ambiguous_read_types_fail(tmp_path):
    r1 = tmp_path / "sample_R1.fastq"
    for second in ("plain", "bad|read_type:unknown"):
        r1.write_text(f"@a|read_type:ont\nAC\n+\nII\n@{second}\nAC\n+\nII\n")
        with pytest.raises(ValueError):
            with assembly_read_libraries(r1):
                pytest.fail("ambiguous types must fail")
