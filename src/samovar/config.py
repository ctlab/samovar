import os
import yaml
import argparse
from typing import Dict, List, Optional, Union
from dataclasses import dataclass
from pathlib import Path

@dataclass
class AnnotatorConfig:
    run_name: str
    type: str
    db_path: str
    db_name: Optional[str] = None
    extra: Optional[str] = None

@dataclass
class PipelineConfig:
    input_config: Optional[str] = None
    input_dir: Optional[str] = None
    output_dir: str = "tests_outs"
    annotators: List[AnnotatorConfig] = None
    read_length: int = 150
    coverage: int = 30
    email: str = "test@samovar.com"

    def __post_init__(self):
        if self.annotators is None:
            self.annotators = []

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'PipelineConfig':
        config = cls()
        
        # Handle input source
        if args.input_config:
            with open(args.input_config, 'r') as f:
                input_config = yaml.safe_load(f)
                config.input_dir = input_config.get('input_dir')
                # Only use config file output_dir if not specified in args
                if not args.output_dir:
                    config.output_dir = input_config.get('output_dir', config.output_dir)
                else:
                    config.output_dir = args.output_dir
                config.read_length = input_config.get('read_length', config.read_length)
                config.coverage = input_config.get('coverage', config.coverage)
                config.email = input_config.get('email', config.email)
                
                # Handle annotators from config
                if 'annotators' in input_config:
                    for ann in input_config['annotators']:
                        config.annotators.append(AnnotatorConfig(**ann))
        elif args.input_dir:
            config.input_dir = args.input_dir
            config.output_dir = args.output_dir

        # Handle command line annotators
        if args.kraken2:
            for k2_config in args.kraken2:
                name, db_path, *extra = k2_config
                config.annotators.append(AnnotatorConfig(
                    run_name=name,
                    type='kraken2',
                    db_path=db_path,
                    extra=' '.join(extra) if extra else None
                ))

        if args.kaiju:
            for kaiju_config in args.kaiju:
                name, db_path, db_name, *extra = kaiju_config
                config.annotators.append(AnnotatorConfig(
                    run_name=name,
                    type='kaiju',
                    db_path=db_path,
                    db_name=db_name,
                    extra=' '.join(extra) if extra else None
                ))

        return config

    def generate_configs(self, base_dir: str) -> Dict[str, str]:
        """Generate all necessary config files and return their paths"""
        base_path = Path(base_dir)
        configs_dir = base_path / '.log' / 'configs'
        configs_dir.mkdir(parents=True, exist_ok=True)
        configs = {}

        # Generate initial annotator config
        init_annotator_config = {
            'r1_dir': str(base_path / 'initial'),
            'r2_dir': str(base_path / 'initial'),
            'output_dir': str(base_path / 'initial_reports'),
            'run_config': [
                {
                    'run_name': ann.run_name,
                    'type': ann.type,
                    'db_path': ann.db_path,
                    **({'db_name': ann.db_name} if ann.db_name else {}),
                    **({'extra': ann.extra} if ann.extra else {})
                }
                for ann in self.annotators
            ]
        }
        init_config_path = configs_dir / 'config_init.yaml'
        with open(init_config_path, 'w') as f:
            yaml.dump(init_annotator_config, f)
        configs['init_annotator'] = str(init_config_path)

        # Generate annotation2iss config
        annotation2iss_config = {
            'annotation_dir': str(base_path / 'initial_annotations'),
            'genome_dir': str(base_path / 'genomes'),
            'output_dir': str(base_path / 'regenerated'),
            'email': self.email,
            'read_length': self.read_length,
            'coverage': self.coverage
        }
        annotation2iss_path = configs_dir / 'config_annotation2iss.yaml'
        with open(annotation2iss_path, 'w') as f:
            yaml.dump(annotation2iss_config, f)
        configs['annotation2iss'] = str(annotation2iss_path)

        # Generate reannotate config
        reannotate_config = {
            'r1_dir': str(base_path / 'regenerated'),
            'r2_dir': str(base_path / 'regenerated'),
            'output_dir': str(base_path / 'regenerated_reports'),
            'run_config': [
                {
                    'run_name': ann.run_name,
                    'type': ann.type,
                    'db_path': ann.db_path,
                    **({'db_name': ann.db_name} if ann.db_name else {}),
                    **({'extra': ann.extra} if ann.extra else {})
                }
                for ann in self.annotators
            ]
        }
        reannotate_path = configs_dir / 'config_reannotate.yaml'
        with open(reannotate_path, 'w') as f:
            yaml.dump(reannotate_config, f)
        configs['reannotate'] = str(reannotate_path)

        return configs

    def generate_pipeline(self, base_dir: str) -> str:
        """Generate the pipeline script and return its path"""
        base_path = Path(base_dir)
        log_dir = base_path / '.log'
        log_dir.mkdir(parents=True, exist_ok=True)
        
        pipeline_path = log_dir / 'samovar.sh'
        
        # Get config paths
        configs = self.generate_configs(base_dir)
        
        # Generate pipeline script
        pipeline_content = f"""# Setup
set -e

if [ -f build/config.json ]; then
    PYTHON_PATH=$(grep -o '"python_path": *"[^"]*"' build/config.json | sed 's/"python_path": *"\\(.*\\)"/\\1/')
    R_PATH=$(grep -o '"r_path": *"[^"]*"' build/config.json | sed 's/"r_path": *"\\(.*\\)"/\\1/')
    R_LIB_PATH=$(grep -o '"r_lib_path": *"[^"]*"' build/config.json | sed 's/"r_lib_path": *"\\(.*\\)"/\\1/')
else
    echo "SamovaR is not installed: check build/config.json"
    exit 1
fi

out_dir="{base_dir}"
mkdir -p $out_dir

# Run annotators on initial reads
snakemake -s workflow/annotators/Snakefile \\
    --configfile {configs['init_annotator']} \\
    --cores 1

# Combine annotation tables
$PYTHON_PATH workflow/combine_annotation_tables.py \\
    -i $out_dir/initial_reports \\
    -o $out_dir/initial_annotations

# Visualize annotations
$R_PATH -s -f "workflow/compare_annotations.R" \\
    --args \\
    --annotation_dir $out_dir/initial_annotations \\
    --output_dir $out_dir/initial_annotations_plots

# Add pre-downloaded genomes to the genome directory
mkdir -p $out_dir/genomes
cp data/test_genomes/meta/* $out_dir/genomes
cp data/test_genomes/host/* $out_dir/genomes

# Translate annotation table to new reads set
snakemake -s workflow/annotation2iss/Snakefile \\
    --configfile {configs['annotation2iss']} \\
    --cores 1

# Clean up
find $out_dir/regenerated -type f -empty -delete
rm $out_dir/regenerated/*_*_*_R*.fastq
rm $out_dir/regenerated/*_abundance*
rm $out_dir/regenerated/*iss.tmp*

# Run annotators on new reads set
snakemake -s workflow/annotators/Snakefile \\
    --configfile {configs['reannotate']} \\
    --cores 1

# Combine annotation tables
$PYTHON_PATH workflow/combine_annotation_tables.py \\
    -i $out_dir/regenerated_reports \\
    -o $out_dir/regenerated_annotations \\
    -s 2

# Visualize & combine results
$R_PATH -s -f "workflow/compare_annotations.R" \\
    --args \\
    --annotation_dir $out_dir/regenerated_annotations \\
    --output_dir $out_dir/regenerated_annotations_plots \\
    --csv $out_dir/regenerated_annotations/combined_annotation_table.csv

# Train and test ML
$PYTHON_PATH workflow/ML.py \\
    --reprofiling_dir $out_dir/initial_annotations \\
    --validation_file $out_dir/regenerated_annotations/combined_annotation_table.csv \\
    --output_dir $out_dir/reprofiled_annotations

# Check reprofiled results
$R_PATH -s -f "workflow/compare_annotations.R" \\
    --args \\
    --annotation_dir $out_dir/reprofiled_annotations \\
    --output_dir $out_dir/reprofiled_annotations_plots \\
    --csv $out_dir/reprofiled_annotations/combined_annotation_table.csv
"""
        
        with open(pipeline_path, 'w') as f:
            f.write(pipeline_content)
        
        return str(pipeline_path)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='SamovaR Pipeline Configuration')
    
    # Input source (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--input-config', help='Path to input configuration file')
    input_group.add_argument('--input-dir', help='Directory containing input FASTQ files')
    
    # Output
    parser.add_argument('--output-dir', default='tests_outs', help='Output directory')
    
    # Annotators
    parser.add_argument('--kraken2', nargs='+', action='append',
                       help='Kraken2 configuration: run_name db_path [extra_args...]')
    parser.add_argument('--kaiju', nargs='+', action='append',
                       help='Kaiju configuration: run_name db_path db_name [extra_args...]')
    
    return parser.parse_args()

def setup_pipeline(args: Optional[argparse.Namespace] = None) -> Dict[str, str]:
    """Main function to set up the pipeline configuration"""
    if args is None:
        args = parse_args()
    
    config = PipelineConfig.from_args(args)
    configs = config.generate_configs(config.output_dir)
    pipeline_path = config.generate_pipeline(config.output_dir)
    
    return {
        'configs': configs,
        'pipeline': pipeline_path
    } 