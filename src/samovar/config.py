import os
import shutil
import sys
import yaml
import argparse
from typing import Dict, List, Optional, Union
from dataclasses import dataclass
from pathlib import Path

from samovar.exec_control import CHECKPOINT_STEPS
from samovar.paths import (
    absolute_path,
    ncbi_email,
    python_path,
    repo_root,
    runtime_path_prefix,
)
from samovar.regenerate import normalize_regeneration_mode

DEFAULT_OUTPUT_DIR = "samovar_out"


def _default_email() -> str:
    return ncbi_email()


def _indent_bash(text: str, spaces: int = 2) -> str:
    pad = " " * spaces
    lines = text.strip("\n").split("\n")
    return "\n".join(pad + line if line.strip() else line for line in lines)


def _checkpoint_block(name: str, body: str) -> str:
    """Wrap a bash snippet so it is skipped when the named checkpoint exists."""
    return (
        f"if ckpt_skip {name}; then\n"
        f"  echo \"[checkpoint] skip {name}\"\n"
        f"else\n"
        f"  echo \"[checkpoint] run {name}\"\n"
        f"{_indent_bash(body, 2)}\n"
        f"  ckpt_finish {name}\n"
        f"fi\n"
    )


def _portable_annotator_cmd(cmd: str) -> str:
    """Keep PATH names in YAML; exec resolves via PATH / install config.

    Absolute binaries the user passed explicitly are stored as-is.
    Warns (does not fail) if the token is not on PATH and not in config.
    """
    from samovar.paths import resolve_executable

    text = (cmd or "").strip()
    if not text:
        return text
    token = text.split()[0]
    path = Path(token).expanduser()
    if path.is_file():
        return text
    resolved = (resolve_executable(token) or "").split()[0]
    if resolved and Path(resolved).expanduser().is_file():
        return text
    if shutil.which(token):
        return text
    if path.is_absolute():
        print(
            f"Warning: annotator binary {token!r} was not found. "
            "Fix the path, or use the command name and set tools / tool_envs "
            "in ~/.config/samovar/config.json so exec can find it.",
            file=sys.stderr,
        )
    else:
        print(
            f"Warning: annotator {token!r} is not on PATH. "
            "Install it, or add its environment in ~/.config/samovar/config.json "
            f"(path, tools.{token}, or tool_envs.{token}) before exec.",
            file=sys.stderr,
        )
    return text


def _absolute_db_path(db_path: str) -> str:
    if not db_path or db_path == ".":
        return db_path or "."
    return absolute_path(db_path)

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
    output_dir: str = DEFAULT_OUTPUT_DIR
    annotators: List[AnnotatorConfig] = None
    read_length: int = 150
    coverage: int = 30
    email: str = None
    cores: int = 1
    max_genomes: int = 50
    regeneration_mode: str = "direct"
    regeneration_n: Optional[int] = None
    regeneration_n_reads: int = 1000
    regeneration_seed: int = 42
    rescale_abundance: bool = False
    gzip_genomes: bool = True
    gzip_reads: bool = False
    reuse_genomes: bool = True
    use_test_genomes: bool = False
    genome_dirs: Optional[List[str]] = None
    run_multiqc: Optional[bool] = None

    def __post_init__(self):
        if self.annotators is None:
            self.annotators = []
        if self.genome_dirs is None:
            self.genome_dirs = []
        if not self.email:
            self.email = _default_email()

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'PipelineConfig':
        config = cls()
        config.cores = getattr(args, 'cores', 1) or 1
        config.max_genomes = getattr(args, 'max_genomes', 50)
        if config.max_genomes is None:
            config.max_genomes = 50
        n_reads = getattr(args, 'N_reads', None)
        if n_reads is None:
            n_reads = getattr(args, 'n_reads', None)
        if n_reads is not None:
            config.regeneration_n_reads = int(n_reads)
        
        # Handle input source
        if args.input_config:
            with open(args.input_config, 'r') as f:
                input_config = yaml.safe_load(f)
                config.input_dir = input_config.get('input_dir')
                cli_out = getattr(args, "output_dir", None)
                config.output_dir = cli_out or input_config.get("output_dir") or DEFAULT_OUTPUT_DIR
                config.read_length = input_config.get('read_length', config.read_length)
                config.coverage = input_config.get('coverage', config.coverage)
                config.email = input_config.get('email', config.email)
                config.regeneration_mode = normalize_regeneration_mode(
                    input_config.get('regeneration_mode', config.regeneration_mode)
                )
                config.regeneration_n = input_config.get('N', config.regeneration_n)
                config.regeneration_n_reads = input_config.get(
                    'N_reads', config.regeneration_n_reads
                )
                config.regeneration_seed = input_config.get('seed', config.regeneration_seed)
                config.rescale_abundance = bool(
                    input_config.get('rescale_abundance', config.rescale_abundance)
                )
                if 'gzip_genomes' in input_config:
                    config.gzip_genomes = bool(input_config.get('gzip_genomes'))
                if 'gzip_reads' in input_config:
                    config.gzip_reads = bool(input_config.get('gzip_reads'))
                
                # Handle annotators from config
                if 'annotators' in input_config:
                    for ann in input_config['annotators']:
                        payload = dict(ann)
                        if payload.get("cmd"):
                            payload["cmd"] = _portable_annotator_cmd(payload["cmd"])
                        if "db_path" in payload:
                            payload["db_path"] = _absolute_db_path(str(payload["db_path"]))
                        config.annotators.append(AnnotatorConfig(**payload))
        elif args.input_dir:
            config.input_dir = args.input_dir
            config.output_dir = getattr(args, "output_dir", None) or DEFAULT_OUTPUT_DIR
        elif getattr(args, "output_dir", None):
            config.output_dir = args.output_dir

        if config.input_dir:
            config.input_dir = absolute_path(config.input_dir)
        config.output_dir = absolute_path(config.output_dir or DEFAULT_OUTPUT_DIR)

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
                        cmd=_portable_annotator_cmd(cmd),
                        db_path=_absolute_db_path(db_path),
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
                        cmd=_portable_annotator_cmd(cmd),
                        db_path=_absolute_db_path(db_path),
                        extra=extra
                    ))

        if getattr(args, "gzip_genomes", None) is not None:
            config.gzip_genomes = bool(args.gzip_genomes)
        if getattr(args, "gzip_reads", None) is not None:
            config.gzip_reads = bool(args.gzip_reads)
        if getattr(args, "reuse_genomes", None) is not None:
            config.reuse_genomes = bool(args.reuse_genomes)
        if getattr(args, "use_test_genomes", False):
            config.use_test_genomes = True
        if getattr(args, "run_multiqc", None) is not None:
            config.run_multiqc = bool(args.run_multiqc)
        extra_dirs = getattr(args, "genome_dirs", None)
        if extra_dirs:
            if isinstance(extra_dirs, str):
                config.genome_dirs = [
                    piece.strip()
                    for piece in extra_dirs.replace(";", ":").split(":")
                    if piece.strip()
                ]
            else:
                config.genome_dirs = [str(p) for p in extra_dirs if str(p).strip()]

        # CLI regeneration depth overrides YAML when both are present.
        n_reads = getattr(args, "N_reads", None)
        if n_reads is None:
            n_reads = getattr(args, "n_reads", None)
        if n_reads is not None:
            config.regeneration_n_reads = int(n_reads)

        return config

    def generate_configs(self, base_dir: str) -> Dict[str, str]:
        """Generate all necessary config files and return their paths"""
        base_path = Path(absolute_path(base_dir))
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
            'regeneration_mode': self.regeneration_mode,
            'N_reads': self.regeneration_n_reads,
            'seed': self.regeneration_seed,
            'rescale_abundance': self.rescale_abundance,
            'gzip_genomes': self.gzip_genomes,
            'gzip_reads': self.gzip_reads,
        }
        if self.regeneration_n:
            annotation2iss_config['N'] = self.regeneration_n
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
        base_dir = absolute_path(base_dir)
        base_path = Path(base_dir)
        log_dir = base_path / '.log'
        log_dir.mkdir(parents=True, exist_ok=True)
        
        pipeline_path = log_dir / 'samovar.sh'
        
        # Get config paths
        configs = self.generate_configs(base_dir)
        root = repo_root()
        py = python_path()
        wf = root / "workflow"
        src = root / "src"
        email = self.email or ncbi_email()
        step_names = " ".join(CHECKPOINT_STEPS)
        tool_path = runtime_path_prefix()
        extra_genome_dirs = ":".join(self.genome_dirs or [])
        reuse_flag = "1" if self.reuse_genomes else "0"
        test_flag = "1" if self.use_test_genomes else "0"
        from samovar.paths import discover_multiqc, shell_source_install_env_snippet

        if self.run_multiqc is None:
            multiqc_flag = "1" if discover_multiqc() else "0"
        else:
            multiqc_flag = "1" if self.run_multiqc else "0"

        env_snippet = shell_source_install_env_snippet()
        
        # Generate pipeline script (absolute paths so exec works from any cwd).
        # Completed steps write $out_dir/.log/checkpoints/<name>.done and are
        # skipped on the next exec unless SAMOVAR_REDO=1 / --redo.
        header = f"""# Setup
set -e
export SAMOVAR_ROOT="{root}"
export NCBI_EMAIL="${{NCBI_EMAIL:-{email}}}"
# Prefer the current install env (other-HPC reinstall) when present.
{env_snippet}# Baked tool bins from config.json: path[], tools.*, tool_envs.*, python/iss.
# Override further with SAMOVAR_PATH=dir1:dir2 (prepended).
export PATH="${{SAMOVAR_PATH:+$SAMOVAR_PATH:}}{tool_path}:{root}/bin:$PATH"
export PYTHONPATH="{src}${{PYTHONPATH:+:$PYTHONPATH}}"
# Honor env PYTHON_PATH (install.sh writes it); fall back to the prepare-time interpreter.
PYTHON_PATH="${{PYTHON_PATH:-{py}}}"
if [ -z "$PYTHON_PATH" ] || [ ! -x "$PYTHON_PATH" ]; then
  PYTHON_PATH="$(command -v python3 || command -v python || true)"
fi
PYTHON_PATH="${{PYTHON_PATH:-python3}}"

out_dir="{base_dir}"
CKPT="$out_dir/.log/checkpoints"
mkdir -p "$CKPT"
# Checkpoint names: {step_names}

ckpt_skip() {{
  [ "${{SAMOVAR_REDO:-0}}" = "0" ] && [ -f "$CKPT/$1.done" ]
}}

ckpt_finish() {{
  "$PYTHON_PATH" -m samovar.exec_control mark "$out_dir" "$1"
  echo "[checkpoint] done $1"
}}

cleanup_tmp_if_requested() {{
  if [ "${{SAMOVAR_CLEANUP_TMP:-0}}" != "0" ]; then
    echo "[cleanup] removing temporary directories under $out_dir"
    "$PYTHON_PATH" -m samovar.exec_control cleanup "$out_dir" || true
  fi
}}

mkdir -p "$out_dir"
mkdir -p "$out_dir/initial" "$out_dir/initial_reports" "$out_dir/regenerated" "$out_dir/regenerated_reports"
# NCBI genome cache reuse (truncated data/test_genomes is not a library).
export SAMOVAR_REUSE_GENOMES="${{SAMOVAR_REUSE_GENOMES:-{reuse_flag}}}"
export SAMOVAR_ALLOW_TEST_GENOMES="${{SAMOVAR_ALLOW_TEST_GENOMES:-{test_flag}}}"
export SAMOVAR_MULTIQC="${{SAMOVAR_MULTIQC:-{multiqc_flag}}}"
export SAMOVAR_GENOME_DIRS="${{SAMOVAR_GENOME_DIRS:+$SAMOVAR_GENOME_DIRS:}}{extra_genome_dirs}"
export SAMOVAR_RUN_DIR="$out_dir"
export SAMOVAR_RUN_GENOMES="$out_dir/genomes"
# All bulky caches live under outdir — never default to ~/.cache (home quota).
export XDG_CACHE_HOME="${{XDG_CACHE_HOME:-$out_dir/.cache}}"
export SAMOVAR_GENOMES="${{SAMOVAR_GENOMES:-$out_dir/.cache/samovar/genomes}}"
export SAMOVAR_PROCESSED_GENOMES="${{SAMOVAR_PROCESSED_GENOMES:-$out_dir/.cache/samovar/genomes}}"
mkdir -p "$XDG_CACHE_HOME/samovar/genomes" "$SAMOVAR_RUN_GENOMES"
"""

        setup_reads = _checkpoint_block(
            "setup_reads",
            f"""# Link/copy source reads into output initial/ when input_dir is provided
# and is not already the destination (generate→preprocess sets input_dir=initial).
if [ -n "{self.input_dir or ''}" ] && [ -d "{self.input_dir or ''}" ]; then
    src_dir=$(readlink -f "{self.input_dir or ''}" || realpath "{self.input_dir or ''}")
    dst_dir=$(readlink -f "$out_dir/initial" || realpath "$out_dir/initial")
    if [ "$src_dir" != "$dst_dir" ]; then
        $PYTHON_PATH -c "from samovar.seqio import link_or_copy_reads; link_or_copy_reads('$src_dir', '$dst_dir')"
    fi
fi""",
        )

        annotate_initial = _checkpoint_block(
            "annotate_initial",
            f"""snakemake -s {wf / 'annotators' / 'Snakefile'} \\
    --configfile {configs['init_annotator']} \\
    --cores {self.cores}""",
        )

        combine_initial = _checkpoint_block(
            "combine_initial",
            f"""$PYTHON_PATH {wf / 'combine_annotation_tables.py'} \\
    -i "$out_dir/initial_reports" \\
    -o "$out_dir/initial_annotations\"""",
        )

        viz_initial = _checkpoint_block(
            "viz_initial",
            f"""$PYTHON_PATH {wf / 'compare_annotations.py'} \\
    --annotation_dir "$out_dir/initial_annotations" \\
    --output_dir "$out_dir/initial_annotations_plots" \\
    --show_top 0""",
        )

        seed_genomes = _checkpoint_block(
            "seed_genomes",
            """# Reuse NCBI/user genome libraries (never the truncated data/test_genomes
# stubs, unless this run's ISS/CAMISIM generate YAML pointed at them or
# SAMOVAR_ALLOW_TEST_GENOMES / --test-genomes is set).
mkdir -p "$out_dir/genomes"
"$PYTHON_PATH" -m samovar.genome_cache seed --dest "$out_dir/genomes" --generate-dir "$out_dir" --genome-dirs "$SAMOVAR_GENOME_DIRS"
""",
        )

        regenerate_reads = _checkpoint_block(
            "regenerate_reads",
            f"""snakemake -s {wf / 'annotation2iss' / 'Snakefile'} \\
    --configfile {configs['annotation2iss']} \\
    --cores {self.cores}
"$PYTHON_PATH" -m samovar.exec_control cleanup "$out_dir" || true""",
        )

        early_exit = """if ! $PYTHON_PATH -c "from samovar.seqio import has_r1_reads; raise SystemExit(0 if has_r1_reads('$out_dir/regenerated') else 1)"; then
    echo "No regenerated reads were produced; skipping re-annotation and reprofiling."
    cleanup_tmp_if_requested
    exit 0
fi
"""

        sort_reads = _checkpoint_block(
            "sort_reads",
            f"""{root / 'bin' / 'samovar'} tools --sort --output_dir "$out_dir\"""",
        )

        annotate_regenerated = _checkpoint_block(
            "annotate_regenerated",
            f"""snakemake -s {wf / 'annotators' / 'Snakefile'} \\
    --configfile {configs['reannotate']} \\
    --cores {self.cores}""",
        )

        combine_regenerated = _checkpoint_block(
            "combine_regenerated",
            f"""$PYTHON_PATH {wf / 'combine_annotation_tables.py'} \\
    -i "$out_dir/regenerated_reports" \\
    -o "$out_dir/regenerated_annotations" \\
    -s 2""",
        )

        viz_regenerated = _checkpoint_block(
            "viz_regenerated",
            f"""$PYTHON_PATH {wf / 'compare_annotations.py'} \\
    --annotation_dir "$out_dir/regenerated_annotations" \\
    --output_dir "$out_dir/regenerated_annotations_plots" \\
    --csv "$out_dir/regenerated_annotations/combined_annotation_table.csv" \\
    --show_top 0""",
        )

        reprofile = _checkpoint_block(
            "reprofile",
            f"""if [ "${{SAMOVAR_ML_FEATURES:-0}}" != "0" ]; then
    echo "[INFO] Extracting per-read features for ML ensemble..."
    $PYTHON_PATH -c "from samovar.seqio import concat_r1_fastqs; concat_r1_fastqs('$out_dir/initial', '$out_dir/combined_temporary_R1.fastq')"
    $PYTHON_PATH {src / 'annotators' / 'fastq_annotator.py'} "$out_dir/combined_temporary_R1.fastq" -o "$out_dir/features.tsv" --chunk_size 50000 --seed {self.regeneration_seed}
    rm -f "$out_dir/combined_temporary_R1.fastq"
    FEATURE_ARG="--features $out_dir/features.tsv"
else
    FEATURE_ARG=""
fi
$PYTHON_PATH {wf / 'ML.py'} \\
    --reprofiling_dir "$out_dir/initial_annotations" \\
    --validation_file "$out_dir/regenerated_annotations/combined_annotation_table.csv" \\
    --output_dir "$out_dir/reprofiled_annotations" \\
    --seed {self.regeneration_seed} \\
    $FEATURE_ARG""",
        )

        viz_reprofiled = _checkpoint_block(
            "viz_reprofiled",
            f"""$PYTHON_PATH {wf / 'compare_annotations.py'} \\
    --annotation_dir "$out_dir/reprofiled_annotations" \\
    --output_dir "$out_dir/reprofiled_annotations_plots" \\
    --csv "$out_dir/reprofiled_annotations/combined_annotation_table.csv" \\
    --show_top 0""",
        )

        footer = """"$PYTHON_PATH" -m samovar.stage_report overview "$out_dir" || true
"$PYTHON_PATH" -m samovar.stage_report bundle "$out_dir" || true
if [ "${SAMOVAR_MULTIQC:-0}" != "0" ]; then
  echo "[multiqc] running MultiQC CLI (--interactive)"
  "$PYTHON_PATH" -m samovar.stage_report multiqc "$out_dir" -- --interactive || true
fi
cleanup_tmp_if_requested
"""

        pipeline_content = "\n".join(
            [
                header,
                setup_reads,
                annotate_initial,
                combine_initial,
                viz_initial,
                seed_genomes,
                regenerate_reads,
                early_exit,
                sort_reads,
                annotate_regenerated,
                combine_regenerated,
                viz_regenerated,
                reprofile,
                viz_reprofiled,
                footer,
            ]
        )

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
    parser.add_argument('--output_dir', default=None, help='Output directory (default: samovar_out, or from --input_config)')
    
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
    try:
        from samovar.genome_cache import remember_prepare_genome_paths

        remember_prepare_genome_paths(extra_dirs=config.genome_dirs)
    except OSError as exc:
        print(f"Warning: could not update genome cache config: {exc}", file=sys.stderr)
    configs = config.generate_configs(config.output_dir)
    pipeline_path = config.generate_pipeline(config.output_dir)
    
    return {
        'configs': configs,
        'pipeline': pipeline_path
    } 