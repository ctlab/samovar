import os
import pytest
import argparse
from pathlib import Path
from src.samovar.iss_config import ISSTestConfig, setup_iss_test
import yaml


@pytest.mark.parametrize("names", [("0", "1"), ("1", "2"), ("patient_A", "batch-02")])
def test_abundance_sample_names_and_count_are_preserved(tmp_path, names):
    from samovar.iss_config import parse_args
    import pandas as pd
    table = tmp_path / "profile.csv"
    pd.DataFrame({"taxid": [562], **{"N_" + name: [1] for name in names}}).to_csv(table, index=False)
    args = parse_args(["--genome_dir", str(tmp_path), "--output_dir", str(tmp_path / "run"), "--abundance", str(table)])
    config = ISSTestConfig.from_args(args)
    assert config.n_samples == 2
    saved = yaml.safe_load(Path(config.generate_config(config.output_dir)).read_text())
    assert saved["sample_names"] == list(names) and saved["n_samples"] == 2
    args.n_samples = 3
    with pytest.raises(ValueError, match="differs"):
        ISSTestConfig.from_args(args)


@pytest.mark.parametrize("header", ["taxid,N_,N_ok", "taxid,N_a/b,N_ok", "taxid,N_a b,N_ok", "taxid,N_a,N_a"])
def test_abundance_rejects_ambiguous_or_unsafe_names_before_generation(tmp_path, header):
    from samovar.abundance import generation_sample_names
    table = tmp_path / "profile.csv"
    table.write_text(header + "\n562,1,1\n")
    with pytest.raises(ValueError, match="Invalid|Duplicate"):
        generation_sample_names(table)

def test_iss_config_from_args(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    
    args = argparse.Namespace(
        genome_dir="/path/to/genomes",
        output_dir=test_output_dir,
        host_genome="/path/to/host.fna",
        n_samples=5,
        total_reads=1000,
        host_fraction="0.1",
        seed=123,
        model="novaseq"
    )
    
    config = ISSTestConfig.from_args(args)
    assert config.genome_dir == str(Path("/path/to/genomes").resolve())
    assert config.output_dir == str(Path(test_output_dir).resolve())
    assert config.host_genome == str(Path("/path/to/host.fna").resolve())
    assert config.n_samples == 5
    assert config.total_reads == 1000
    assert config.host_fraction == "0.1"
    assert config.seed == 123
    assert config.model == "novaseq"

def test_iss_config_defaults(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    
    args = argparse.Namespace(
        genome_dir="/path/to/genomes",
        output_dir=test_output_dir,
        host_genome="/path/to/host.fna",
        n_samples=None,
        total_reads=None,
        host_fraction=None,
        seed=None,
        model=None
    )
    
    config = ISSTestConfig.from_args(args)
    assert config.genome_dir == str(Path("/path/to/genomes").resolve())
    assert config.output_dir == str(Path(test_output_dir).resolve())
    assert config.host_genome == str(Path("/path/to/host.fna").resolve())
    assert config.n_samples == 10  # default
    assert config.total_reads == 2000  # default
    assert config.host_fraction == "RANDOM"  # default
    assert config.seed == 42  # default
    assert config.model == "hiseq"  # default

def test_generate_config(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    
    config = ISSTestConfig(
        genome_dir="/path/to/genomes",
        output_dir=test_output_dir,
        host_genome="/path/to/host.fna",
        n_samples=5,
        total_reads=1000,
        host_fraction="0.1",
        seed=123,
        model="novaseq"
    )
    
    config_path = config.generate_config(test_output_dir)
    assert os.path.exists(config_path)
    
    with open(config_path, 'r') as f:
        config_content = yaml.safe_load(f)
        assert config_content['genome_dir'] == "/path/to/genomes"
        assert config_content['output_dir'] == str((Path(test_output_dir) / 'initial').resolve())
        assert config_content['host_genome'] == "/path/to/host.fna"
        assert config_content['n_samples'] == 5
        assert config_content['total_reads'] == 1000
        assert config_content['host_fraction'] == "0.1"
        assert config_content['seed'] == 123
        assert config_content['model'] == "novaseq"
        assert config_content['genomes'] == []
        assert config_content['max_genomes'] == 100.0

def test_generate_pipeline(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    
    config = ISSTestConfig(
        genome_dir="/path/to/genomes",
        output_dir=test_output_dir,
        host_genome="/path/to/host.fna"
    )
    
    pipeline_path = config.generate_pipeline(test_output_dir)
    assert os.path.exists(pipeline_path)
    
    with open(pipeline_path, 'r') as f:
        pipeline_content = f.read()
        assert "set -e" in pipeline_content
        assert "_runtime-env" in pipeline_content
        assert f"out_dir=\"{Path(test_output_dir).resolve()}\"" in pipeline_content
        assert '"$SNAKEMAKE" -s ' in pipeline_content
        assert '--directory "$out_dir"' in pipeline_content
        assert "workflow/iss_test/Snakefile" in pipeline_content
        assert "samovar/env" in pipeline_content
        assert "_runtime-env" in pipeline_content

def test_setup_iss_test(tmp_path):
    test_output_dir = str(tmp_path / "run")
    os.makedirs(test_output_dir, exist_ok=True)
    
    args = argparse.Namespace(
        genome_dir="/path/to/genomes",
        output_dir=test_output_dir,
        host_genome="/path/to/host.fna",
        n_samples=5,
        total_reads=1000,
        host_fraction="0.1",
        seed=123,
        model="novaseq"
    )
    
    result = setup_iss_test(args)
    assert 'config' in result
    assert 'pipeline' in result
    
    config_path = result['config']
    pipeline_path = result['pipeline']
    
    assert os.path.exists(config_path)
    assert os.path.exists(pipeline_path)
    
    # Verify config content
    with open(config_path, 'r') as f:
        config_content = yaml.safe_load(f)
        assert config_content['genome_dir'] == "/path/to/genomes"
        assert config_content['output_dir'] == str((Path(test_output_dir) / 'initial').resolve())
        assert config_content['host_genome'] == "/path/to/host.fna"
        assert config_content['n_samples'] == 5
        assert config_content['total_reads'] == 1000
        assert config_content['host_fraction'] == "0.1"
        assert config_content['seed'] == 123
        assert config_content['model'] == "novaseq"
    
    # Verify pipeline content
    with open(pipeline_path, 'r') as f:
        pipeline_content = f.read()
        assert "set -e" in pipeline_content
        assert "_runtime-env" in pipeline_content
        assert f"out_dir=\"{Path(test_output_dir).resolve()}\"" in pipeline_content
        assert '"$SNAKEMAKE" -s ' in pipeline_content
        assert '--directory "$out_dir"' in pipeline_content
        assert "workflow/iss_test/Snakefile" in pipeline_content


def test_parse_args_directory_alias():
    from samovar.iss_config import parse_args

    args = parse_args(["--directory", "out", "--genome_dir", "g"])
    assert args.output_dir == "out"


def test_parse_args_reindex_and_raw_genomes():
    from samovar.iss_config import parse_args

    args = parse_args(
        [
            "--output_dir",
            "out",
            "--accessions",
            "GCF_000819615.1",
            "GCF_000840245.1",
            "--reindex",
            "1",
            "--raw-genomes",
            "0",
        ]
    )
    assert args.accessions == ["GCF_000819615.1", "GCF_000840245.1"]
    assert args.reindex == 1
    assert args.raw_genomes == 0


def test_parse_args_reindex_defaults_to_database():
    from samovar.iss_config import parse_args

    args = parse_args(["--output_dir", "out"])
    assert args.reindex == 1


def test_parse_args_genome_skip_and_max_mb():
    from samovar.iss_config import parse_args

    args = parse_args(
        [
            "--output_dir",
            "out",
            "--genome_dir",
            "genomes",
            "--max-genome-mb",
            "inf",
            "--genome-skip-list",
            "9606,562",
        ]
    )
    assert args.max_genome_mb == float("inf")
    assert args.genome_skip_list == "9606,562"


def test_iss_config_writes_skip_list(tmp_path):
    from samovar.iss_config import ISSTestConfig

    cfg = ISSTestConfig(
        genome_dir="/path/to/genomes",
        output_dir=str(tmp_path),
        host_genome="",
        max_genome_mb=12.5,
        genome_skip_list="9606",
    )
    path = cfg.generate_config(str(tmp_path))
    data = yaml.safe_load(Path(path).read_text())
    assert data["max_genome_mb"] == 12.5
    assert data["genome_skip_list"] == "9606"


def test_iss_config_max_genomes_cli(tmp_path):
    from samovar.iss_config import parse_args

    args = parse_args(
        [
            "--output_dir",
            str(tmp_path),
            "--genome_dir",
            str(tmp_path),
            "--max-genomes",
            "2",
        ]
    )
    assert args.max_genomes == 2.0
    cfg = ISSTestConfig.from_args(args)
    assert cfg.max_genomes == 2.0
    data = yaml.safe_load(Path(cfg.generate_config(str(tmp_path))).read_text())
    assert data["max_genomes"] == 2.0
    inf_args = parse_args(
        ["--output_dir", str(tmp_path), "--max-genomes", "inf"]
    )
    assert inf_args.max_genomes == float("inf")


def test_parse_args_threads_maps_to_iss_cpus():
    from samovar.iss_config import parse_args

    args = parse_args(
        [
            "--output_dir",
            "out",
            "--threads",
            "8",
        ]
    )
    assert args.cores == 8
    assert "- - c" not in (args.extra_flags or "")

@pytest.mark.parametrize("field", ["reads_generator", "metagenome_generator"])
def test_generate_preflight_runs_before_custom_tool_validation(tmp_path, monkeypatch, field):
    from samovar.iss_config import parse_args
    import samovar.prepare_import as preflight
    import samovar.reads_generators as reads
    import samovar.metagenome_generators as meta

    args = parse_args(["--output_dir", str(tmp_path / "out"),
                       "--genome_dir", str(tmp_path), "--" + field, "custom_new"])
    calls = []
    monkeypatch.setattr(preflight, "offer_prepare_imports", lambda ns: calls.append(("offer", getattr(ns, field))))

    def validate(name):
        calls.append(("validate", name))
        raise ValueError("stop after preflight")

    module = reads if field == "reads_generator" else meta
    monkeypatch.setattr(module, "require_known_" + field, validate)
    with pytest.raises(ValueError, match="stop after preflight"):
        setup_iss_test(args)
    assert calls == [("offer", "custom_new"), ("validate", "custom_new")]


def test_generate_autoimport_flag():
    from samovar.iss_config import parse_args

    assert parse_args(["--output_dir", "out", "--autoimport", "0"]).autoimport == 0
    assert parse_args(["--output_dir", "out", "--autoimport", "1"]).autoimport == 1
