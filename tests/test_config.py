import os
import pytest
import argparse
from pathlib import Path
from src.samovar.config import PipelineConfig, setup_pipeline, AnnotatorConfig
import yaml

def test_pipeline_config_from_args_file(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    config_path = Path(test_output_dir) / "test_config.yaml"
    config_content = """
    input_dir: /path/to/input
    output_dir: /path/to/output
    read_length: 100
    coverage: 20
    email: test@example.com
    annotators:
      - run_name: k2-test
        type: kraken2
        cmd: /path/to/kraken2
        db_path: /path/to/kraken_db
        extra: --minK 2 --maxK 10
      - run_name: kaiju-test
        type: kaiju
        cmd: /path/to/kaiju
        db_path: /path/to/kaiju
    """
    with open(config_path, 'w') as f:
        f.write(config_content)
    args = argparse.Namespace(
        input_config=str(config_path),
        input_dir=None,
        output_dir=test_output_dir,
        kraken2=None,
        kaiju=None
    )
    config = PipelineConfig.from_args(args)
    assert config.input_dir == str(Path("/path/to/input").resolve())
    assert config.output_dir == str(Path(test_output_dir).resolve())
    assert config.read_length == 100
    assert config.coverage == 20
    assert config.email == "test@example.com"
    assert len(config.annotators) == 2
    k2 = next(a for a in config.annotators if a.type == "kraken2")
    assert k2.run_name == "k2-test"
    assert k2.cmd == "/path/to/kraken2"
    assert k2.db_path == "/path/to/kraken_db"
    assert k2.extra == "--minK 2 --maxK 10"
    kaiju = next(a for a in config.annotators if a.type == "kaiju")
    assert kaiju.run_name == "kaiju-test"
    assert kaiju.cmd == "/path/to/kaiju"
    assert kaiju.db_path == "/path/to/kaiju"

def test_pipeline_config_uses_yaml_output_dir_when_cli_unset(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    config_path = Path(test_output_dir) / "test_config.yaml"
    config_path.write_text(
        """
input_dir: /path/to/input
output_dir: /from/yaml/out
annotators: []
"""
    )
    args = argparse.Namespace(
        input_config=str(config_path),
        input_dir=None,
        output_dir=None,
        kraken2=None,
        kaiju=None,
    )
    config = PipelineConfig.from_args(args)
    assert config.output_dir == str(Path("/from/yaml/out").resolve())


def test_pipeline_config_default_output_dir():
    args = argparse.Namespace(
        input_config=None,
        input_dir="/path/to/input",
        output_dir=None,
        kraken2=None,
        kaiju=None,
    )
    config = PipelineConfig.from_args(args)
    assert config.output_dir == str((Path.cwd() / "samovar_out").resolve())

def test_pipeline_config_from_args_cli(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    args = argparse.Namespace(
        input_config=None,
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        kraken2=[["/path/to/kraken2 /path/to/kraken_db --minK 2 --maxK 10"]],
        kaiju=[["/path/to/kaiju /path/to/kaiju"]]
    )
    config = PipelineConfig.from_args(args)
    assert config.input_dir == str(Path("/path/to/input").resolve())
    assert config.output_dir == str(Path(test_output_dir).resolve())
    assert len(config.annotators) == 2
    k2 = next(a for a in config.annotators if a.type == "kraken2")
    assert k2.run_name == "kraken2"
    assert k2.cmd == "/path/to/kraken2"
    assert k2.db_path == "/path/to/kraken_db"
    assert k2.extra == "--minK 2 --maxK 10"
    kaiju = next(a for a in config.annotators if a.type == "kaiju")
    assert kaiju.run_name == "kaiju"
    assert kaiju.cmd == "/path/to/kaiju"
    assert kaiju.db_path == "/path/to/kaiju"

def test_prepare_rejects_colliding_sample_names(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "sample_R1.fastq").write_text("@a\nA\n+\nI\n")
    (reads / "sample_1.fq").write_text("@b\nT\n+\nI\n")
    with pytest.raises(ValueError, match="Multiple R1 FASTQ files"):
        PipelineConfig.from_args(argparse.Namespace(
            input_config=None, input_dir=str(reads), output_dir=str(tmp_path / "out"),
        ))
    assert not (tmp_path / "out").exists()


def test_prepare_rejects_collision_with_existing_run(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    first, second, out = (tmp_path / name for name in ("first", "second", "out"))
    first.mkdir()
    second.mkdir()
    (first / "sample_R1.fastq").write_text("@a\nA\n+\nI\n")
    (second / "sample_R1.fastq").write_text("@b\nT\n+\nI\n")
    (out / "initial").mkdir(parents=True)
    (out / "initial" / "sample_R1.fastq").symlink_to(first / "sample_R1.fastq")
    with pytest.raises(ValueError, match="Sample name collision"):
        PipelineConfig.from_args(argparse.Namespace(
            input_config=None, input_dir=str(second), output_dir=str(out),
        ))
    assert not (out / ".log").exists()


def test_generate_configs(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    config = PipelineConfig(
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        annotators=[
            AnnotatorConfig(
                run_name="k2-test",
                type="kraken2",
                cmd="/path/to/kraken2",
                db_path="/path/to/kraken_db",
                extra="--minK 2 --maxK 10"
            ),
            AnnotatorConfig(
                run_name="kaiju-test",
                type="kaiju",
                cmd="/path/to/kaiju",
                db_path="/path/to/kaiju",
                db_name="kaiju_db.fmi"
            )
        ]
    )
    configs = config.generate_configs(test_output_dir)
    configs_dir = Path(test_output_dir) / '.log' / 'configs'
    assert configs_dir.exists()
    assert configs_dir.is_dir()
    assert os.path.exists(configs['init_annotator'])
    assert os.path.exists(configs['annotation2iss'])
    assert os.path.exists(configs['reannotate'])
    assert os.path.exists(configs['scoring'])
    with open(configs['init_annotator'], 'r') as f:
        init_config = yaml.safe_load(f)
        assert init_config['r1_dir'] == str((Path(test_output_dir) / 'initial_trimmed').resolve())
        assert init_config['r2_dir'] == str((Path(test_output_dir) / 'initial_trimmed').resolve())
        assert init_config['output_dir'] == str((Path(test_output_dir) / 'initial_reports').resolve())
        assert len(init_config['run_config']) == 2
        k2_config = next(c for c in init_config['run_config'] if c['type'] == 'kraken2')
        assert k2_config['cmd'] == "/path/to/kraken2"
        kaiju_config = next(c for c in init_config['run_config'] if c['type'] == 'kaiju')
        assert kaiju_config['cmd'] == "/path/to/kaiju"
    with open(configs['annotation2iss'], 'r') as f:
        iss_config = yaml.safe_load(f)
        assert iss_config['gzip_genomes'] is True
        assert iss_config['gzip_reads'] is False
        assert iss_config['N_reads'] == 1000000
        assert iss_config['max_genomes'] == 100
        assert iss_config['max_genome_mb'] == 100
        assert iss_config['genome_skip_list'] == []


def test_pipeline_config_n_reads_cli_overrides_default(tmp_path):
    test_output_dir = str(tmp_path / "run")
    args = argparse.Namespace(
        input_config=None,
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        kraken2=None,
        kaiju=None,
        N_reads=20000,
        max_genomes=20,
    )
    config = PipelineConfig.from_args(args)
    assert config.regeneration_n_reads == 20000
    assert config.max_genomes == 20
    configs = config.generate_configs(test_output_dir)
    with open(configs["annotation2iss"], "r") as f:
        iss_config = yaml.safe_load(f)
    assert iss_config["N_reads"] == 20000
    assert iss_config["max_genomes"] == 20


def test_pipeline_config_max_genomes_default(tmp_path, monkeypatch):
    monkeypatch.delenv("SAMOVAR_MAX_GENOMES", raising=False)
    args = argparse.Namespace(
        input_config=None,
        input_dir="/path/to/input",
        output_dir=str(tmp_path),
        kraken2=None,
        kaiju=None,
    )
    config = PipelineConfig.from_args(args)
    assert config.max_genomes == 100
    configs = config.generate_configs(str(tmp_path))
    with open(configs["annotation2iss"], "r") as f:
        iss_config = yaml.safe_load(f)
    assert iss_config["max_genomes"] == 100


def test_pipeline_config_genome_skip_and_max_mb_cli(tmp_path):
    args = argparse.Namespace(
        input_config=None,
        input_dir="/path/to/input",
        output_dir=str(tmp_path),
        kraken2=None,
        kaiju=None,
        max_genome_mb=50,
        genome_skip_list="9606,562",
    )
    config = PipelineConfig.from_args(args)
    assert config.max_genome_mb == 50
    assert "9606" in config.genome_skip_list
    assert "562" in config.genome_skip_list
    configs = config.generate_configs(str(tmp_path))
    with open(configs["annotation2iss"], "r") as f:
        iss_config = yaml.safe_load(f)
    assert iss_config["max_genome_mb"] == 50
    assert "9606" in iss_config["genome_skip_list"]


def test_generate_pipeline(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    config = PipelineConfig(
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        annotators=[
            AnnotatorConfig(
                run_name="k2-test",
                type="kraken2",
                cmd="/path/to/kraken2",
                db_path="/path/to/kraken_db",
                extra="--minK 2 --maxK 10"
            ),
            AnnotatorConfig(
                run_name="kaiju-test",
                type="kaiju",
                cmd="/path/to/kaiju",
                db_path="/path/to/kaiju",
                db_name="kaiju_db.fmi"
            )
        ]
    )
    
    # Test pipeline generation
    pipeline_path = config.generate_pipeline(test_output_dir)
    assert os.path.exists(pipeline_path)
    
    # Read and verify pipeline content
    with open(pipeline_path, 'r') as f:
        pipeline_content = f.read()
        
        # Check basic structure
        assert "set -e" in pipeline_content
        assert "_runtime-env" in pipeline_content
        assert "workflow/compare_annotations.py" in pipeline_content
        assert "workflow/compare_annotations.R" not in pipeline_content
        assert "export R_PATH=" not in pipeline_content
        
        # Check config file paths
        assert f"out_dir=\"{Path(test_output_dir).resolve()}\"" in pipeline_content
        assert "_runtime-env" in pipeline_content
        assert '[ -f "$CKPT/$1.done" ]' in pipeline_content
        assert "visualization failed; continuing" not in pipeline_content
        assert "XDG_CONFIG_HOME" in pipeline_content
        assert "samovar/env" in pipeline_content
        assert "_runtime-env" in pipeline_content
        assert 'export PATH="$_SAMOVAR_CALLER_PATH"' in pipeline_content
        
        # Check snakemake commands
        assert '"$SNAKEMAKE" -s ' in pipeline_content
        assert '--directory "$out_dir"' in pipeline_content
        assert "workflow/annotators/Snakefile" in pipeline_content
        assert "workflow/abundance2iss/Snakefile" in pipeline_content
        
        # Check Python and R commands
        assert "workflow/combine_annotation_tables.py" in pipeline_content
        assert "samovar.combine_tables finalize" in pipeline_content
        init_yaml = Path(test_output_dir) / ".log" / "configs" / "config_init.yaml"
        assert "incremental_annotations: true" in init_yaml.read_text(encoding="utf-8")
        assert "workflow/compare_annotations.py" in pipeline_content
        assert "workflow/ML.py" in pipeline_content
        assert "--seed 42" in pipeline_content
        assert "samovar.seqio" in pipeline_content
        assert "link_or_copy_reads" in pipeline_content
        assert "validate_fastq_directory" in pipeline_content
        assert "cp data/test_genomes" not in pipeline_content
        assert "samovar.genome_cache seed" in pipeline_content
        assert "export PYTHONPATH=" in pipeline_content
        assert "export SAMOVAR_ROOT=" in pipeline_content
        assert "samovar.exec_control" in pipeline_content
        assert "ckpt_skip" in pipeline_content
        assert "ckpt_finish" in pipeline_content
        assert 'SAMOVAR_START="${SAMOVAR_START:-setup_reads}"' in pipeline_content
        assert 'SAMOVAR_END="${SAMOVAR_END:-viz_reprofiled}"' in pipeline_content
        assert "exec_control require" in pipeline_content
        assert "exec_control active" in pipeline_content
        assert "samovar.stage_report" in pipeline_content
        assert "samovar.scorers run" in pipeline_content
        assert "SAMOVAR_MULTIQC" in pipeline_content
        for step in (
            "setup_reads",
            "qc_initial",
            "annotate_initial",
            "combine_initial",
            "viz_initial",
            "abundance_tables",
            "regenerate_tables",
            "score_regenerated_tables",
            "seed_genomes",
            "regenerate_reads",
            "sort_reads",
            "qc_generated",
            "annotate_regenerated",
            "combine_regenerated",
            "viz_regenerated",
            "reprofile",
            "viz_reprofiled",
        ):
            assert f"ckpt_skip {step}" in pipeline_content
            assert f"ckpt_finish {step}" in pipeline_content


def test_generate_pipeline_multiqc_flags(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    annotators = [
        AnnotatorConfig(
            run_name="k2-test",
            type="kraken2",
            cmd="/path/to/kraken2",
            db_path="/path/to/kraken_db",
        )
    ]
    off = PipelineConfig(
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        annotators=annotators,
        run_multiqc=False,
    )
    off_text = Path(off.generate_pipeline(test_output_dir)).read_text()
    assert 'export SAMOVAR_MULTIQC="${SAMOVAR_MULTIQC:-0}"' in off_text
    on = PipelineConfig(
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        annotators=annotators,
        run_multiqc=True,
    )
    on_text = Path(on.generate_pipeline(test_output_dir)).read_text()
    assert 'export SAMOVAR_MULTIQC="${SAMOVAR_MULTIQC:-1}"' in on_text
    assert "stage_report multiqc" in on_text


def test_setup_pipeline(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    
    args = argparse.Namespace(
        input_config=None,
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        kraken2=[["/path/to/kraken2 /path/to/kraken_db --minK 2 --maxK 10"]],
        kaiju=[["/path/to/kaiju /path/to/kaiju kaiju_db.fmi"]]
    )

    result = setup_pipeline(args)
    
    # Check that all expected files exist
    assert 'configs' in result
    assert 'pipeline' in result
    
    configs = result['configs']
    pipeline_path = result['pipeline']
    
    assert os.path.exists(configs['init_annotator'])
    assert os.path.exists(configs['annotation2iss'])
    assert os.path.exists(configs['reannotate'])
    assert os.path.exists(pipeline_path)
    
    # Verify pipeline content
    with open(pipeline_path, 'r') as f:
        pipeline_content = f.read()
        assert f"out_dir=\"{Path(test_output_dir).resolve()}\"" in pipeline_content
        assert '--directory "$out_dir"' in pipeline_content


def test_generate_pipeline_resolves_tools_without_baking_env_bins(tmp_path, monkeypatch):
    env = tmp_path / "kaiju_env"
    (env / "bin").mkdir(parents=True)
    kaiju = env / "bin" / "kaiju"
    kaiju.write_text("#!/bin/sh\n")
    kaiju.chmod(0o755)
    extra = tmp_path / "other_env" / "bin"
    extra.mkdir(parents=True)
    monkeypatch.setattr(
        "samovar.paths.load_config",
        lambda: {
            "path": [str(extra.parent)],
            "tool_envs": {"kaiju": str(env)},
            "tools": {"kaiju": str(kaiju)},
        },
    )
    out = tmp_path / "out"
    (tmp_path / "reads").mkdir()
    config = PipelineConfig(
        input_dir=str(tmp_path / "reads"),
        output_dir=str(out),
        annotators=[
            AnnotatorConfig(
                run_name="kaiju-test",
                type="kaiju",
                cmd="kaiju",
                db_path=str(tmp_path / "db"),
            )
        ],
    )
    script = Path(config.generate_pipeline(str(out))).read_text()
    assert str(env / "bin") not in script
    assert str(extra) not in script
    assert "samovar/env" in script
    assert "_runtime-env" in script


def test_prepare_to_flags_embed_convert(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    out = tmp_path / "out"
    args = argparse.Namespace(
        input_config=None,
        input_dir=str(reads),
        output_dir=str(out),
        kraken2=[["kraken2 /tmp/k2"]],
        kaiju=None,
        to_mpa=True,
        to_cami=True,
        export_to=["kraken2"],
        export_taxdump="/tmp/taxdump",
        export_taxonomy="gtdb",
    )
    result = setup_pipeline(args)
    text = Path(result["pipeline"]).read_text()
    assert "samovar.abundance_correctors" in text
    assert "--to cami" in text
    assert "--to kraken2_mpa" in text
    assert "--to kraken2" in text
    assert "exports/initial" in text
    assert "exports/regenerated" in text
    assert "exports/reprofiled" in text
    assert "--taxdump /tmp/taxdump" in text
    assert "--taxonomy gtdb" in text
    cfg = PipelineConfig.from_args(args)
    assert cfg.export_formats == ["kraken2", "kraken2_mpa", "cami"]
    assert cfg.export_corrector == "logistic"


def test_prepare_without_to_flags_defaults_logistic_export(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    args = argparse.Namespace(
        input_config=None,
        input_dir=str(reads),
        output_dir=str(tmp_path / "out"),
        kraken2=[["kraken2 /tmp/k2"]],
        kaiju=None,
    )
    text = Path(setup_pipeline(args)["pipeline"]).read_text()
    assert "samovar.abundance_correctors" in text
    assert "--to abundance" in text
    assert "--export logistic" in text
    cfg = PipelineConfig.from_args(args)
    assert cfg.export_corrector == "logistic"
    assert cfg.export_formats == ["abundance"]


def test_prepare_no_export_skips_corrector(tmp_path):
    reads = tmp_path / "reads"
    reads.mkdir()
    args = argparse.Namespace(
        input_config=None,
        input_dir=str(reads),
        output_dir=str(tmp_path / "out"),
        kraken2=[["kraken2 /tmp/k2"]],
        kaiju=None,
        skip_export=True,
    )
    text = Path(setup_pipeline(args)["pipeline"]).read_text()
    assert "samovar.abundance_correctors" not in text
    cfg = PipelineConfig.from_args(args)
    assert cfg.export_corrector == "off"


def test_prepare_yaml_export_formats(tmp_path):
    cfg = tmp_path / "p.yaml"
    cfg.write_text(
        "input_dir: /in\nexport_formats: [cami, mpa]\ntaxonomy: gtdb\nannotators: []\n",
        encoding="utf-8",
    )
    args = argparse.Namespace(
        input_config=str(cfg),
        input_dir=None,
        output_dir=str(tmp_path / "out"),
        kraken2=None,
        kaiju=None,
    )
    config = PipelineConfig.from_args(args)
    assert config.export_formats == ["cami", "kraken2_mpa"]
    assert config.export_taxonomy == "gtdb"


def test_prepare_preserves_tagged_hybrid_generator(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "s_R1.fastq").write_text(
        "@short|taxid:562|read_type:illumina\nACGT\n+\nIIII\n"
        "@long|taxid:562|read_type:ont\nACGTACGT\n+\nIIIIIIII\n"
    )
    args = argparse.Namespace(input_config=None, input_dir=str(reads), output_dir=str(tmp_path / "out"))
    assert PipelineConfig.from_args(args).reads_generator == "hybrid"
    args.reads_generator = "iss"
    assert PipelineConfig.from_args(args).reads_generator == "iss"
    args.reads_generator = None
    config = tmp_path / "input.yaml"
    config.write_text(yaml.safe_dump({"input_dir": str(reads)}))
    args.input_config = str(config)
    assert PipelineConfig.from_args(args).reads_generator == "hybrid"
    config.write_text(yaml.safe_dump({"input_dir": str(reads), "reads_generator": "iss"}))
    assert PipelineConfig.from_args(args).reads_generator == "iss"


@pytest.mark.parametrize("mode, expected", [("hybrid", "hybrid"), ("ont", "nanosim"), ("illumina", "art"), ("wgsim", "wgsim"), ("table", "camisim")])
def test_prepare_inherits_pending_hybrid_generation(tmp_path, monkeypatch, mode, expected):
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    out = tmp_path / "run"
    generated = out / ".generate" / "configs" / "camisim.yaml"
    generated.parent.mkdir(parents=True)
    generated.write_text(f"mode: {mode}\n")
    args = argparse.Namespace(input_config=None, input_dir=str(out / "initial"), output_dir=str(out))
    assert PipelineConfig.from_args(args).reads_generator == expected
    args.reads_generator = "iss"
    assert PipelineConfig.from_args(args).reads_generator == "iss"


@pytest.mark.parametrize("field,generator", [("reads_generator", g) for g in ["iss", "nanosim", "art", "wgsim", "hybrid"]] + [("metagenome_generator", "nanosim"), ("metagenome_generator", "hybrid")])
def test_prepare_inherits_generate_snapshot(tmp_path, monkeypatch, field, generator):
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    source = tmp_path / "source"
    snapshot = source / ".hydra" / "config.yaml"
    snapshot.parent.mkdir(parents=True)
    snapshot.write_text(yaml.safe_dump({"generate": {"args": {field: generator}}}))
    # A stale CAMISIM YAML must not override the last recorded generate command.
    stale = source / ".generate" / "configs" / "camisim.yaml"
    stale.parent.mkdir(parents=True)
    stale.write_text("mode: hybrid\n")
    args = argparse.Namespace(input_config=None, input_dir=str(source / "initial"), output_dir=str(tmp_path / "new"))
    assert getattr(PipelineConfig.from_args(args), field) == generator
    args.reads_generator = "iss"
    assert PipelineConfig.from_args(args).reads_generator == "iss"


@pytest.mark.parametrize("field", ["reads_generator", "metagenome_generator"])
@pytest.mark.parametrize("legacy", [False, True])
def test_prepare_inherits_camisim_ont_mode(tmp_path, monkeypatch, field, legacy):
    monkeypatch.setenv("SAMOVAR_CONFIG", str(tmp_path / "config.json"))
    if legacy:
        saved = tmp_path / ".generate" / "configs" / "camisim.yaml"
        payload = {field: "camisim", "mode": "ont"}
    else:
        saved = tmp_path / ".hydra" / "config.yaml"
        payload = {"generate": {"args": {field: "camisim", "camisim_mode": "ont"}}}
    saved.parent.mkdir(parents=True)
    saved.write_text(yaml.safe_dump(payload))
    args = argparse.Namespace(input_config=None, input_dir=str(tmp_path / "initial"), output_dir=str(tmp_path))
    assert getattr(PipelineConfig.from_args(args), field) == "nanosim"
