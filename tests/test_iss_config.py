import os
import shutil
import pytest
import argparse
from pathlib import Path
from src.samovar.iss_config import ISSTestConfig, setup_iss_test
import yaml

def clean_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)

def test_iss_config_from_args():
    test_output_dir = 'tests_outs/test_iss_config_from_args'
    clean_dir(test_output_dir)
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

def test_iss_config_defaults():
    test_output_dir = 'tests_outs/test_iss_config_defaults'
    clean_dir(test_output_dir)
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

def test_generate_config():
    test_output_dir = 'tests_outs/test_generate_config'
    clean_dir(test_output_dir)
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
        assert config_content['max_genomes'] == float("inf")

def test_generate_pipeline():
    test_output_dir = 'tests_outs/test_generate_pipeline'
    clean_dir(test_output_dir)
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
        assert "PYTHON_PATH=" in pipeline_content
        assert f"out_dir=\"{Path(test_output_dir).resolve()}\"" in pipeline_content
        assert "snakemake -s " in pipeline_content
        assert "workflow/iss_test/Snakefile" in pipeline_content
        assert "samovar/env" in pipeline_content
        assert "SAMOVAR_PATH" in pipeline_content

def test_setup_iss_test():
    test_output_dir = 'tests_outs/test_setup_iss_test'
    clean_dir(test_output_dir)
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
        assert "PYTHON_PATH=" in pipeline_content
        assert f"out_dir=\"{Path(test_output_dir).resolve()}\"" in pipeline_content
        assert "snakemake -s " in pipeline_content
        assert "workflow/iss_test/Snakefile" in pipeline_content


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