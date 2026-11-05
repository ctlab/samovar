from pathlib import Path

import pytest

from samovar.iss_config import parse_args, setup_iss_test
from samovar.parse_annotators import extract_true_taxid
from samovar.reads_generators import (
    MissingReadsGeneratorError,
    extra_ids_for_generator,
    flags_apply_to_reads_generator,
    get_reads_generator,
    reads_output_stems,
    require_known_reads_generator,
    resolve_reads_generator,
)
from samovar.seqio import list_fastq_samples
from samovar.tools_import import import_tool


def test_resolve_builtin_reads_generators():
    assert resolve_reads_generator("iss") == ("builtin", "iss")
    assert resolve_reads_generator("art_illumina") == ("builtin", "art")
    assert resolve_reads_generator("wgsim") == ("builtin", "wgsim")
    assert resolve_reads_generator("nanosim3") == ("builtin", "nanosim")
    assert resolve_reads_generator("hybrid") == ("builtin", "hybrid")
    assert resolve_reads_generator("myreads") == ("custom", "myreads")


def test_output_stems_hybrid_and_prepare():
    assert reads_output_stems("1", extra_ids=["illumina"]) == ["1_illumina", "1_full"]
    assert reads_output_stems(
        "1", annotator="kaiju", extra_ids=["sequence_type"], stage="regenerate"
    ) == ["1_kaiju_sequence_type", "1_kaiju"]


def test_list_fastq_samples_drops_tech_when_full_present(tmp_path):
    (tmp_path / "1_full_R1.fastq").write_text("@a\nA\n+\nI\n")
    (tmp_path / "1_illumina_R1.fastq").write_text("@b\nA\n+\nI\n")
    (tmp_path / "1_kaiju_R1.fastq").write_text("@c\nA\n+\nI\n")
    (tmp_path / "1_kaiju_ont_R1.fastq").write_text("@d\nA\n+\nI\n")
    samples = list_fastq_samples(tmp_path)
    assert "1_full" in samples
    assert "1_illumina" not in samples
    assert "1_kaiju" in samples
    assert "1_kaiju_ont" not in samples


def test_missing_custom_reads_generator(tmp_path, monkeypatch):
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    with pytest.raises(MissingReadsGeneratorError, match="ghost_reads"):
        require_known_reads_generator("ghost_reads")


def test_custom_python_generate_and_flags(tmp_path, monkeypatch):
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    script = tmp_path / "myreads.py"
    script.write_text(
        "from pathlib import Path\n"
        "from samovar.camisim import tag_read_id\n"
        "\n"
        "def generate(spec, metadata, config):\n"
        "    out = Path(spec['output_dir'])\n"
        "    out.mkdir(parents=True, exist_ok=True)\n"
        "    sample = '1'\n"
        "    extra = list(spec.get('extra_ids') or [])\n"
        "    paths = []\n"
        "    for stem in ([f'{sample}_{e}' for e in extra] + [f'{sample}_full']):\n"
        "        r1 = out / f'{stem}_R1.fastq'\n"
        "        r2 = out / f'{stem}_R2.fastq'\n"
        "        hid = tag_read_id('read0', extra[0] if extra else '', '562')\n"
        "        rec = f'@{hid}\\nACGT\\n+\\nIIII\\n'\n"
        "        r1.write_text(rec)\n"
        "        r2.write_text(rec)\n"
        "        paths.extend([str(r1), str(r2)])\n"
        "        (out / 'argv.txt').write_text(' '.join(config.get('extra_argv') or []))\n"
        "    return paths\n"
    )
    import_tool(
        name="myreads",
        tool_type="reads",
        exec_path=str(script),
        flags="--profile hiseq",
        also_repo_build=False,
    )
    out = tmp_path / "initial"
    gen = get_reads_generator("myreads")
    written = gen.generate(
        None,
        {
            "output_dir": str(out),
            "genome_dir": str(tmp_path),
            "host_fraction": "0",
            "n_samples": 1,
            "total_reads": 10,
            "reads_generator_flags": "--keep-tmp",
            "filename_ids": "sequence_type",
            "stage": "generate",
        },
    )
    assert any(Path(p).name == "1_full_R1.fastq" for p in written)
    assert any(Path(p).name == "1_sequence_type_R1.fastq" for p in written)
    text = (out / "1_full_R1.fastq").read_text()
    assert extract_true_taxid(text.split()[0]) == "562"
    assert "read_type:sequence_type" in text
    argv = (out / "argv.txt").read_text()
    assert "--profile" in argv
    assert "hiseq" in argv
    assert "--keep-tmp" in argv


def test_prepare_flags_reads_generator(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    binary = tmp_path / "myreads.py"
    binary.write_text("def generate(spec, metadata, config):\n    return []\n")
    import_tool(
        name="myreads",
        tool_type="reads",
        exec_path=str(binary),
        flags="--from-import",
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
            "reads_generator": "myreads",
            "tool_flags": [["reads_generator", "--from-prepare"]],
        },
    )()
    config = PipelineConfig.from_args(args)
    assert config.reads_generator == "myreads"
    assert "--from-prepare" in (config.reads_generator_flags or "")
    text = Path(config.generate_configs(str(tmp_path / "out"))["annotation2iss"]).read_text()
    assert "myreads" in text
    assert "reads_generator_flags" in text


def test_prepare_unknown_reads_generator(tmp_path, monkeypatch):
    from samovar.config import PipelineConfig
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    (tmp_path / "reads").mkdir()
    args = type(
        "Args",
        (),
        {
            "input_config": None,
            "input_dir": str(tmp_path / "reads"),
            "output_dir": str(tmp_path / "out"),
            "reads_generator": "ghost_reads",
        },
    )()
    with pytest.raises(MissingReadsGeneratorError, match="ghost_reads"):
        PipelineConfig.from_args(args)


def test_generate_cli_flags_and_abundance(tmp_path):
    args = parse_args(
        [
            "--genome_dir",
            str(tmp_path),
            "--output_dir",
            str(tmp_path / "run"),
            "--host_genome",
            "",
            "--reads_generator",
            "iss",
            "-i",
            str(tmp_path / "abund.csv"),
            "--flags",
            "iss",
            "--gc_bias",
        ]
    )
    assert args.reads_generator == "iss"
    assert args.abundance_table.endswith("abund.csv")
    assert args.tool_flags == [["iss", "--gc_bias"]]


def test_setup_generate_merges_iss_flags(tmp_path, monkeypatch):
    from samovar.paths import write_config

    cfg = tmp_path / "config.json"
    monkeypatch.setenv("SAMOVAR_CONFIG", str(cfg))
    write_config({"root": str(tmp_path), "tools": {}}, also_repo_build=False)
    genomes = tmp_path / "genomes"
    genomes.mkdir()
    (genomes / "562.fna").write_text(">g\nACGT\n")
    args = parse_args(
        [
            "--genome_dir",
            str(genomes),
            "--output_dir",
            str(tmp_path / "run"),
            "--flags",
            "reads_generator",
            "--gc_bias",
        ]
    )
    result = setup_iss_test(args)
    data = Path(result["config"]).read_text()
    assert "gc_bias" in data
    assert "iss" in data


def test_flags_apply_reads():
    assert flags_apply_to_reads_generator("iss", "iss")
    assert flags_apply_to_reads_generator("reads_generator", "art")
    assert flags_apply_to_reads_generator("art", "illumina")
    assert flags_apply_to_reads_generator("nanosim", "ont")
    assert extra_ids_for_generator("art") == ["illumina"]
    assert extra_ids_for_generator("wgsim") == ["wgsim"]
    assert extra_ids_for_generator("nanosim") == ["ont"]


@pytest.mark.parametrize("generator", ["nanosim", "ont", "art", "wgsim", "hybrid", "meta:nanosim"])
@pytest.mark.parametrize("gzip_reads", [False, True])
def test_regeneration_uses_selected_simulator(tmp_path, monkeypatch, generator, gzip_reads):
    import pandas as pd

    from samovar import camisim, table2iss
    from samovar.reads_generators import simulate_from_sample_tables
    from samovar.seqio import fastq_pair_paths, iter_fastq_records

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    monkeypatch.setattr(camisim, "discover_camisim", lambda: None)
    monkeypatch.setattr(camisim, "camisim_taxdump", lambda root: "")
    monkeypatch.setattr(camisim, "available_camisim_types", lambda root: ["art", "nanosim3", "wgsim"])

    def unexpected_iss(**kwargs):
        pytest.fail("Non-ISS regeneration must not launch ISS")

    monkeypatch.setattr(table2iss, "generate_reads_metagenome", unexpected_iss)
    monkeypatch.setattr(camisim, "design_communities", unexpected_iss)
    calls = []

    def simulate(cfg, work, cam_type, dest):
        communities = camisim.read_distribution_files(work / "distributions")
        calls.append((cam_type, communities, cfg["extra_flags"]))
        dest.mkdir(parents=True)
        length = 4000 if cam_type == "nanosim3" else 150
        text = "".join(
            f"@{taxid}_read\n{'A' * length}\n+\n{'I' * length}\n"
            for taxid, amount in communities[0].items() if amount > 0
        )
        (dest / "sample_0_01.fastq").write_text(text)
        if cam_type != "nanosim3":
            (dest / "sample_0_02.fastq").write_text(text)
        return True

    monkeypatch.setattr(camisim, "_run_camisim_reads", simulate)
    genomes = {}
    for taxid in ("562", "1280"):
        # The resolved taxid must survive paths whose basename is an accession.
        genome = tmp_path / f"GCF_{taxid}.fna"
        genome.write_text(">contig\n" + "ACGT" * 1000 + "\n")
        genomes[taxid] = str(genome)
    tables = {
        "sample_a": pd.DataFrame({"taxid": ["562", "1280"], "N_kaiju": [9, 1]}),
        "sample_b": pd.DataFrame({"taxid": ["562", "1280"], "N_kaiju": [1, 9]}),
    }
    out = tmp_path / "regenerated"
    config = {"reads_generator_flags": "--custom-option value"}
    if generator.startswith("meta:"):
        config = {"metagenome_generator": generator.split(":")[1],
                  "metagenome_generator_flags": "--custom-option value"}
    if generator == "hybrid":
        out.mkdir()
        for sample in tables:
            stale, _ = fastq_pair_paths(out / f"{sample}_kaiju", gzip_reads=not gzip_reads)
            from samovar.seqio import open_text
            with open_text(stale, "wt") as handle:
                handle.write("@stale|taxid:562|read_type:illumina\nACGT\n+\nIIII\n")
    simulate_from_sample_tables(
        tables, genomes, str(out), annotator_cols={"kaiju": "N_kaiju"},
        reads_generator="iss" if generator.startswith("meta:") else generator,
        total_amount=20, gzip_reads=gzip_reads, config=config,
    )
    expected_types = {
        "nanosim": ["nanosim3"], "ont": ["nanosim3"], "art": ["art"],
        "wgsim": ["wgsim"], "hybrid": ["art", "nanosim3"], "meta:nanosim": ["nanosim3"],
    }[generator]
    assert [call[0] for call in calls] == expected_types * 2
    assert calls[0][1] == [{"562": 0.9, "1280": 0.1}]
    assert calls[-1][1] == [{"562": 0.1, "1280": 0.9}]
    assert all(call[2] == "--custom-option value" for call in calls)
    for sample in tables:
        r1, r2 = fastq_pair_paths(out / f"{sample}_kaiju", gzip_reads=gzip_reads)
        if generator == "hybrid":
            from samovar.seqio import list_fastq_samples
            assert sample + "_kaiju" in list_fastq_samples(out)
        records = list(iter_fastq_records(r1))
        lengths = {len(rec[1].strip()) for rec in records}
        assert lengths == ({150, 4000} if generator == "hybrid" else
                           {4000} if generator in ("nanosim", "ont", "meta:nanosim") else {150})
        assert {extract_true_taxid(rec[0]) for rec in records} == {"562", "1280"}
        assert Path(r2).exists() == (generator in ("art", "wgsim"))


def test_regeneration_simulator_failure_does_not_fall_back_to_iss(tmp_path, monkeypatch):
    import pandas as pd

    from samovar import camisim
    from samovar.reads_generators import simulate_from_sample_tables

    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    monkeypatch.setattr(camisim, "default_generate_config", lambda **kwargs: {})

    def fail(config):
        raise RuntimeError("NanoSim failed")

    monkeypatch.setattr(camisim, "run_from_config", fail)
    with pytest.raises(RuntimeError, match="NanoSim failed"):
        simulate_from_sample_tables(
            {"sample": pd.DataFrame({"taxid": ["562"], "N_kaiju": [10]})},
            {"562": "genome.fna"}, str(tmp_path / "regenerated"),
            annotator_cols={"kaiju": "N_kaiju"}, reads_generator="nanosim",
        )
