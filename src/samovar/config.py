import os
import yaml
import argparse
from typing import Dict, List, Optional, Union
from dataclasses import dataclass
from pathlib import Path

def _default_email() -> str:
    for key in ("NCBI_EMAIL", "ENTREZ_EMAIL", "SAMOVAR_EMAIL"):
        value = os.environ.get(key, "").strip()
        if value:
            return value
    return "anonymous@example.com"

@dataclass
class AnnotatorConfig:
    run_name: str
    type: str
    db_path: str
    cmd: str
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
    email: str = None
    cores: int = 1
    max_genomes: int = 50

    def __post_init__(self):
        if self.annotators is None:
            self.annotators = []
        if not self.email:
            self.email = _default_email()

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'PipelineConfig':
        config = cls()
        config.cores = getattr(args, 'cores', 1) or 1
        config.max_genomes = getattr(args, 'max_genomes', 50)
        if config.max_genomes is None:
            config.max_genomes = 50
        
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
        # Get all attributes that start with 'cmd_'
        cmd_attrs = [attr for attr in dir(args) if attr.startswith('cmd_')]
        
        for attr in cmd_attrs:
            if getattr(args, attr):
                for cmd_config in getattr(args, attr):
                    # Split the command string into parts (db is optional for dummy/custom tools)
                    parts = cmd_config[0].split()
                    if len(parts) < 1:
                        raise ValueError(f"Invalid command format for {attr}: {cmd_config[0]}")
                    
                    # First part is the command path
                    cmd = parts[0]
                    # Second part is the database path (optional)
                    db_path = parts[1] if len(parts) > 1 else "."
                    # Any remaining parts are extra arguments
                    extra = ' '.join(parts[2:]) if len(parts) > 2 else None
                    
                    # Extract type from command basename
                    cmd_basename = os.path.basename(cmd)
                    type_name = cmd_basename.split('.')[0]
                    
                    # Extract run name from attribute name (remove 'cmd_' prefix)
                    run_name = attr[4:]
                    dummy_aliases = {"dummy", "dummy9606", "constant9606", "constant", "random"}
                    if run_name.lower() in dummy_aliases or type_name.lower() in dummy_aliases:
                        type_name = "constant9606"
                    # Konstantaza CLI: --custom-test "metauto /path/to/db"
                    if run_name.lower() in {"custom", "custom-test", "custom_test"}:
                        run_name = type_name
                    
                    config.annotators.append(AnnotatorConfig(
                        run_name=run_name,
                        type=type_name,
                        cmd=cmd,
                        db_path=db_path,
                        extra=extra
                    ))

        # Handle legacy command line annotators
        for attr in dir(args):
            if attr in ['kraken2', 'kaiju', 'dummy'] and getattr(args, attr):
                for cmd_config in getattr(args, attr):
                    # Split the command string into parts (db is optional for dummy)
                    parts = cmd_config[0].split()
                    if len(parts) < 1:
                        raise ValueError(f"Invalid command format for {attr}: {cmd_config[0]}")
                    
                    # First part is the command path
                    cmd = parts[0]
                    # Second part is the database path (optional)
                    db_path = parts[1] if len(parts) > 1 else "."
                    # Any remaining parts are extra arguments
                    extra = ' '.join(parts[2:]) if len(parts) > 2 else None
                    
                    # Extract type from command basename
                    cmd_basename = os.path.basename(cmd)
                    type_name = cmd_basename.split('.')[0]
                    dummy_aliases = {"dummy", "dummy9606", "constant9606", "constant", "random"}
                    if attr == "dummy" or type_name.lower() in dummy_aliases:
                        type_name = "constant9606"
                    
                    config.annotators.append(AnnotatorConfig(
                        run_name=type_name if attr != "dummy" else "dummy",
                        type=type_name,
                        cmd=cmd,
                        db_path=db_path,
                        extra=extra
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
                    'cmd': ann.cmd,
                    'db_path': ann.db_path,
                    'threads': self.cores,
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
            'coverage': self.coverage,
            'max_genomes': getattr(self, 'max_genomes', 50),
            'cores': self.cores,
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
                    'cmd': ann.cmd,
                    'db_path': ann.db_path,
                    'threads': self.cores,
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
mkdir -p $out_dir/initial $out_dir/initial_reports $out_dir/regenerated $out_dir/regenerated_reports

# Link/copy source reads into output initial/ when input_dir is provided
# and is not already the destination (generate→preprocess sets input_dir=initial).
if [ -n "{self.input_dir}" ] && [ -d "{self.input_dir}" ]; then
    src_dir=$(readlink -f "{self.input_dir}")
    dst_dir=$(readlink -f "$out_dir/initial")
    if [ "$src_dir" != "$dst_dir" ]; then
        ln -sf $(readlink -f {self.input_dir}/*_R*.fastq) $out_dir/initial/ 2>/dev/null || cp {self.input_dir}/*_R*.fastq $out_dir/initial/
    fi
fi

# Run annotators on initial reads
snakemake -s workflow/annotators/Snakefile \\
    --configfile {configs['init_annotator']} \\
    --cores {self.cores}

# Combine annotation tables
$PYTHON_PATH workflow/combine_annotation_tables.py \\
    -i $out_dir/initial_reports \\
    -o $out_dir/initial_annotations

# Visualize annotations
$R_PATH -s -f "workflow/compare_annotations.R" \\
    --args \\
    --annotation_dir $out_dir/initial_annotations \\
    --output_dir $out_dir/initial_annotations_plots \\
    --show_top 0

# Add pre-downloaded genomes to the genome directory
mkdir -p $out_dir/genomes
cp data/test_genomes/meta/* $out_dir/genomes
cp data/test_genomes/host/* $out_dir/genomes

# Translate annotation table to new reads set
snakemake -s workflow/annotation2iss/Snakefile \\
    --configfile {configs['annotation2iss']} \\
    --cores {self.cores}

# Clean up
{{
    find $out_dir/regenerated -type f -empty -delete || true
    rm -f $out_dir/regenerated/*processed* || true
    rm -f $out_dir/regenerated/*_abundance* || true
    rm -f $out_dir/regenerated/*iss.tmp* || true
}} || {{
    echo "Warning: Some cleanup operations failed"
}}

if ! ls $out_dir/regenerated/*_R1.fastq >/dev/null 2>&1; then
    echo "No regenerated reads were produced; skipping re-annotation and reprofiling."
    exit 0
fi

# Sort paired-end reads to ensure matching order
bin/samovar tools --sort --output_dir $out_dir

# Run annotators on new reads set
snakemake -s workflow/annotators/Snakefile \\
    --configfile {configs['reannotate']} \\
    --cores {self.cores}

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
    --csv $out_dir/regenerated_annotations/combined_annotation_table.csv \\
    --show_top 0

# Train and test ML
if [ "${{SAMOVAR_ML_FEATURES:-0}}" != "0" ]; then
    echo "[INFO] Extracting per-read features for ML ensemble..."
    cat $out_dir/initial/*_R1.fastq > $out_dir/combined_temporary_R1.fastq
    $PYTHON_PATH src/annotators/fastq_annotator.py $out_dir/combined_temporary_R1.fastq -o $out_dir/features.tsv --chunk_size 50000
    rm -f $out_dir/combined_temporary_R1.fastq
    FEATURE_ARG="--features $out_dir/features.tsv"
else
    FEATURE_ARG=""
fi
$PYTHON_PATH workflow/ML.py \\
    --reprofiling_dir $out_dir/initial_annotations \\
    --validation_file $out_dir/regenerated_annotations/combined_annotation_table.csv \\
    --output_dir $out_dir/reprofiled_annotations \\
    $FEATURE_ARG

# Check reprofiled results
$R_PATH -s -f "workflow/compare_annotations.R" \\
    --args \\
    --annotation_dir $out_dir/reprofiled_annotations \\
    --output_dir $out_dir/reprofiled_annotations_plots \\
    --csv $out_dir/reprofiled_annotations/combined_annotation_table.csv \\
    --show_top 0
"""
        
        with open(pipeline_path, 'w') as f:
            f.write(pipeline_content)
        
        return str(pipeline_path)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='SamovaR Pipeline Configuration')
    
    # Input source (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=False)
    input_group.add_argument('--input_config', help='Path to input configuration file', required=False)
    input_group.add_argument('--input_dir', help='Directory containing input FASTQ files', required=False)
    
    # Output
    parser.add_argument('--output_dir', default='tests_outs', help='Output directory', required=True)
    
    # Add a dynamic command argument group
    cmd_group = parser.add_argument_group('Command Arguments')
    cmd_group.add_argument('--cmd_*', action='append', nargs=1, metavar='CMD',
                          help='Command specification in format: "cmd_path db_path [extra_args...]"')
    
    args = parser.parse_args()
    
    # Process command arguments
    cmd_args = {}
    for arg, value in vars(args).items():
        if arg.startswith('cmd_'):
            cmd_args[arg] = value
    
    # Remove the original cmd_* argument
    if hasattr(args, 'cmd_*'):
        delattr(args, 'cmd_*')
    
    # Add processed command arguments
    for arg, value in cmd_args.items():
        setattr(args, arg, value)
    
    return args

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