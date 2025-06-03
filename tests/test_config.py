import os
import shutil
import pytest
import argparse
from pathlib import Path
from src.samovar.config import PipelineConfig, setup_pipeline, AnnotatorConfig
import yaml

def clean_dir(path):
    if os.path.exists(path):
        shutil.rmtree(path)

def test_pipeline_config_from_args_file():
    test_output_dir = 'tests_outs/test_pipeline_config_from_args_file'
    clean_dir(test_output_dir)
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
        db_path: /path/to/kraken_db
        extra: --minK 2 --maxK 10
      - run_name: kaiju-test
        type: kaiju
        db_path: /path/to/kaiju
        db_name: kaiju_db.fmi
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
    assert config.input_dir == "/path/to/input"
    assert config.output_dir == test_output_dir
    assert config.read_length == 100
    assert config.coverage == 20
    assert config.email == "test@example.com"
    assert len(config.annotators) == 2
    k2 = next(a for a in config.annotators if a.type == "kraken2")
    assert k2.run_name == "k2-test"
    assert k2.db_path == "/path/to/kraken_db"
    assert k2.extra == "--minK 2 --maxK 10"
    kaiju = next(a for a in config.annotators if a.type == "kaiju")
    assert kaiju.run_name == "kaiju-test"
    assert kaiju.db_path == "/path/to/kaiju"
    assert kaiju.db_name == "kaiju_db.fmi"

def test_pipeline_config_from_args_cli():
    test_output_dir = 'tests_outs/test_pipeline_config_from_args_cli'
    clean_dir(test_output_dir)
    os.makedirs(test_output_dir, exist_ok=True)
    args = argparse.Namespace(
        input_config=None,
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        kraken2=[["k2-test", "/path/to/kraken_db", "--minK", "2", "--maxK", "10"]],
        kaiju=[["kaiju-test", "/path/to/kaiju", "kaiju_db.fmi"]]
    )
    config = PipelineConfig.from_args(args)
    assert config.input_dir == "/path/to/input"
    assert config.output_dir == test_output_dir
    assert len(config.annotators) == 2
    k2 = next(a for a in config.annotators if a.type == "kraken2")
    assert k2.run_name == "k2-test"
    assert k2.db_path == "/path/to/kraken_db"
    assert k2.extra == "--minK 2 --maxK 10"
    kaiju = next(a for a in config.annotators if a.type == "kaiju")
    assert kaiju.run_name == "kaiju-test"
    assert kaiju.db_path == "/path/to/kaiju"
    assert kaiju.db_name == "kaiju_db.fmi"

def test_generate_configs():
    test_output_dir = 'tests_outs/test_generate_configs'
    clean_dir(test_output_dir)
    os.makedirs(test_output_dir, exist_ok=True)
    config = PipelineConfig(
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        annotators=[
            AnnotatorConfig(
                run_name="k2-test",
                type="kraken2",
                db_path="/path/to/kraken_db",
                extra="--minK 2 --maxK 10"
            ),
            AnnotatorConfig(
                run_name="kaiju-test",
                type="kaiju",
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
    with open(configs['init_annotator'], 'r') as f:
        init_config = yaml.safe_load(f)
        assert init_config['r1_dir'] == str(Path(test_output_dir) / 'initial')
        assert init_config['r2_dir'] == str(Path(test_output_dir) / 'initial')
        assert init_config['output_dir'] == str(Path(test_output_dir) / 'initial_reports')
        assert len(init_config['run_config']) == 2
    with open(configs['annotation2iss'], 'r') as f:
        iss_config = yaml.safe_load(f)
        assert iss_config['annotation_dir'] == str(Path(test_output_dir) / 'initial_annotations')
        assert iss_config['genome_dir'] == str(Path(test_output_dir) / 'genomes')
        assert iss_config['output_dir'] == str(Path(test_output_dir) / 'regenerated')
        assert iss_config['read_length'] == 150
        assert iss_config['coverage'] == 30
    with open(configs['reannotate'], 'r') as f:
        reannotate_config = yaml.safe_load(f)
        assert reannotate_config['r1_dir'] == str(Path(test_output_dir) / 'regenerated')
        assert reannotate_config['r2_dir'] == str(Path(test_output_dir) / 'regenerated')
        assert reannotate_config['output_dir'] == str(Path(test_output_dir) / 'regenerated_reports')
        assert len(reannotate_config['run_config']) == 2

def test_generate_pipeline():
    test_output_dir = 'tests_outs/test_generate_pipeline'
    clean_dir(test_output_dir)
    os.makedirs(test_output_dir, exist_ok=True)
    config = PipelineConfig(
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        annotators=[
            AnnotatorConfig(
                run_name="k2-test",
                type="kraken2",
                db_path="/path/to/kraken_db",
                extra="--minK 2 --maxK 10"
            ),
            AnnotatorConfig(
                run_name="kaiju-test",
                type="kaiju",
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
        assert "PYTHON_PATH=" in pipeline_content
        assert "R_PATH=" in pipeline_content
        
        # Check config file paths
        assert f"out_dir=\"{test_output_dir}\"" in pipeline_content
        
        # Check snakemake commands
        assert "snakemake -s workflow/annotators/Snakefile" in pipeline_content
        assert "snakemake -s workflow/annotation2iss/Snakefile" in pipeline_content
        
        # Check Python and R commands
        assert "workflow/combine_annotation_tables.py" in pipeline_content
        assert "workflow/compare_annotations.R" in pipeline_content
        assert "workflow/ML.py" in pipeline_content

def test_setup_pipeline():
    test_output_dir = 'tests_outs/test_setup_pipeline'
    clean_dir(test_output_dir)
    os.makedirs(test_output_dir, exist_ok=True)
    
    args = argparse.Namespace(
        input_config=None,
        input_dir="/path/to/input",
        output_dir=test_output_dir,
        kraken2=[["k2-test", "/path/to/kraken_db", "--minK", "2", "--maxK", "10"]],
        kaiju=[["kaiju-test", "/path/to/kaiju", "kaiju_db.fmi"]]
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
        assert f"out_dir=\"{test_output_dir}\"" in pipeline_content 