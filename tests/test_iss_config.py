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
    assert config.genome_dir == "/path/to/genomes"
    assert config.output_dir == test_output_dir
    assert config.host_genome == "/path/to/host.fna"
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
    assert config.genome_dir == "/path/to/genomes"
    assert config.output_dir == test_output_dir
    assert config.host_genome == "/path/to/host.fna"
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
        assert config_content['output_dir'] == str(Path(test_output_dir) / 'initial')
        assert config_content['host_genome'] == "/path/to/host.fna"
        assert config_content['n_samples'] == 5
        assert config_content['total_reads'] == 1000
        assert config_content['host_fraction'] == "0.1"
        assert config_content['seed'] == 123
        assert config_content['model'] == "novaseq"
        assert config_content['genomes'] == []

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
        assert f"out_dir=\"{test_output_dir}\"" in pipeline_content
        assert "snakemake -s workflow/iss_test/Snakefile" in pipeline_content

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
        assert config_content['output_dir'] == str(Path(test_output_dir) / 'initial')
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
        assert f"out_dir=\"{test_output_dir}\"" in pipeline_content
        assert "snakemake -s workflow/iss_test/Snakefile" in pipeline_content 