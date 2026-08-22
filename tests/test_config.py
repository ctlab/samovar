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

def test_pipeline_config_uses_yaml_output_dir_when_cli_unset():
    test_output_dir = "tests_outs/test_pipeline_config_yaml_output_dir"
    clean_dir(test_output_dir)
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

def test_pipeline_config_from_args_cli():
    test_output_dir = 'tests_outs/test_pipeline_config_from_args_cli'
    clean_dir(test_output_dir)
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
    with open(configs['init_annotator'], 'r') as f:
        init_config = yaml.safe_load(f)
        assert init_config['r1_dir'] == str((Path(test_output_dir) / 'initial').resolve())
        assert init_config['r2_dir'] == str((Path(test_output_dir) / 'initial').resolve())
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
        assert "PYTHON_PATH=" in pipeline_content
        assert "workflow/compare_annotations.py" in pipeline_content
        assert "workflow/compare_annotations.R" not in pipeline_content
        assert "R_PATH=" not in pipeline_content
        
        # Check config file paths
        assert f"out_dir=\"{Path(test_output_dir).resolve()}\"" in pipeline_content
        assert 'PYTHON_PATH="${PYTHON_PATH:-' in pipeline_content
        assert '[ -f "$CKPT/$1.done" ]' in pipeline_content
        assert "visualization failed; continuing" not in pipeline_content
        
        # Check snakemake commands
        assert "snakemake -s " in pipeline_content
        assert "workflow/annotators/Snakefile" in pipeline_content
        assert "workflow/annotation2iss/Snakefile" in pipeline_content
        
        # Check Python and R commands
        assert "workflow/combine_annotation_tables.py" in pipeline_content
        assert "workflow/compare_annotations.py" in pipeline_content
        assert "workflow/ML.py" in pipeline_content
        assert "samovar.seqio" in pipeline_content
        assert "link_or_copy_reads" in pipeline_content
        assert "cp data/test_genomes" not in pipeline_content
        assert "export PYTHONPATH=" in pipeline_content
        assert "export SAMOVAR_ROOT=" in pipeline_content
        assert "samovar.exec_control" in pipeline_content
        assert "ckpt_skip" in pipeline_content
        assert "ckpt_finish" in pipeline_content
        for step in (
            "setup_reads",
            "annotate_initial",
            "combine_initial",
            "viz_initial",
            "seed_genomes",
            "regenerate_reads",
            "sort_reads",
            "annotate_regenerated",
            "combine_regenerated",
            "viz_regenerated",
            "reprofile",
            "viz_reprofiled",
        ):
            assert f"ckpt_skip {step}" in pipeline_content
            assert f"ckpt_finish {step}" in pipeline_content

def test_setup_pipeline():
    test_output_dir = 'tests_outs/test_setup_pipeline'
    clean_dir(test_output_dir)
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