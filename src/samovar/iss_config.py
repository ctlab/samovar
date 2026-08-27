import os
import yaml
import argparse
from typing import Dict, Optional
from dataclasses import dataclass
from pathlib import Path

@dataclass
class ISSTestConfig:
    genome_dir: str
    output_dir: str
    host_genome: str
    n_samples: int = 10
    total_reads: int = 2000
    host_fraction: str = "RANDOM"
    seed: int = 42
    model: str = "hiseq"
    cores: int = 1
    extra_flags: str = ""
    reads_generator: str = "iss"
    abundance_table: str = ""

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'ISSTestConfig':
        from samovar.paths import absolute_path

        abundance = getattr(args, "abundance_table", None) or ""
        if abundance:
            abundance = absolute_path(abundance)
        extra = getattr(args, "extra_flags", None) or ""
        return cls(
            genome_dir=absolute_path(args.genome_dir) if args.genome_dir else "",
            output_dir=absolute_path(args.output_dir),
            host_genome=absolute_path(args.host_genome) if args.host_genome else "",
            n_samples=args.n_samples if args.n_samples is not None else 10,
            total_reads=args.total_reads if args.total_reads is not None else 2000,
            host_fraction=args.host_fraction if args.host_fraction is not None else "RANDOM",
            seed=args.seed if args.seed is not None else 42,
            model=args.model if args.model is not None else "hiseq",
            cores=args.cores if getattr(args, 'cores', None) is not None else 1,
            extra_flags=str(extra),
            reads_generator=str(getattr(args, "reads_generator", None) or "iss"),
            abundance_table=str(abundance or ""),
        )

    def generate_config(self, base_dir: str) -> str:
        """Generate ISS test config file and return its path"""
        from samovar.paths import absolute_path

        base_path = Path(absolute_path(base_dir))
        configs_dir = base_path / '.generate' / 'configs'
        configs_dir.mkdir(parents=True, exist_ok=True)
        
        config_path = configs_dir / 'iss_config.yaml'
        
        config = {
            'genome_dir': self.genome_dir,
            'output_dir': str(base_path / 'initial'),
            'host_genome': self.host_genome,
            'n_samples': self.n_samples,
            'total_reads': self.total_reads,
            'host_fraction': self.host_fraction,
            'seed': self.seed,
            'model': self.model,
            'cores': self.cores,
            'genomes': [],  # Will be automatically populated
            'reads_generator': self.reads_generator,
        }
        if self.extra_flags:
            config['extra_flags'] = self.extra_flags
            config['reads_generator_flags'] = self.extra_flags
        if self.abundance_table:
            config['abundance_table'] = self.abundance_table
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        return str(config_path)

    def generate_pipeline(self, base_dir: str) -> str:
        """Generate the ISS test pipeline script and return its path"""
        from samovar.paths import absolute_path, python_path, repo_root, runtime_path_prefix

        base_dir = absolute_path(base_dir)
        base_path = Path(base_dir)
        generate_dir = base_path / '.generate'
        generate_dir.mkdir(parents=True, exist_ok=True)
        
        pipeline_path = generate_dir / 'generate.sh'
        
        # Get config path
        config_path = self.generate_config(base_dir)
        
        # Generate pipeline script
        root = repo_root()
        py = python_path()
        tool_path = runtime_path_prefix()
        snakefile = root / "workflow" / "iss_test" / "Snakefile"
        from samovar.paths import shell_source_install_env_snippet

        env_snippet = shell_source_install_env_snippet()
        pipeline_content = f"""# Setup
set -e
export SAMOVAR_ROOT="{root}"
{env_snippet}export PATH="${{SAMOVAR_PATH:+$SAMOVAR_PATH:}}{tool_path}:{root}/bin:$PATH"
export PYTHONPATH="{root / 'src'}${{PYTHONPATH:+:$PYTHONPATH}}"
PYTHON_PATH="${{PYTHON_PATH:-{py}}}"
if [ -z "$PYTHON_PATH" ] || [ ! -x "$PYTHON_PATH" ]; then
  PYTHON_PATH="$(command -v python3 || command -v python || true)"
fi
PYTHON_PATH="${{PYTHON_PATH:-python3}}"

out_dir="{base_dir}"
mkdir -p "$out_dir"

# Generate reads with InSilicoSeq
snakemake -s {snakefile} \\
    --configfile {config_path} \\
    --cores {self.cores}
"""
        
        with open(pipeline_path, 'w') as f:
            f.write(pipeline_content)
        
        return str(pipeline_path)

def parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='SamovaR ISS Test Configuration')
    
    # Required arguments
    parser.add_argument('--genome_dir', required=False, default=None,
                       help='Directory containing reference genomes')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for generated files')
    parser.add_argument('--host_genome', required=False, default="",
                       help='Path to host genome file')
    parser.add_argument(
        '--accessions',
        nargs='+',
        default=None,
        help='NCBI assembly accessions (GCF_*/GCA_*) to download for this generate',
    )
    parser.add_argument(
        '--reindex',
        type=int,
        default=0,
        choices=(0, 1, 2),
        help='0: download to $out/.genomes (no index); 1: samovar_database + index; '
             '2: $out/.genomes + index',
    )
    parser.add_argument(
        '--raw-genomes',
        type=int,
        default=0,
        choices=(0, 1),
        help='Keep raw NCBI FASTA after parsing (default 0: delete raw)',
    )
    
    # Optional arguments with defaults
    parser.add_argument('--n_samples', type=int, default=10,
                       help='Number of samples to generate')
    parser.add_argument('--total_reads', type=int, default=2000,
                       help=(
                           'Sequencing records per sample: paired fragments for '
                           'ISS/Illumina/wgsim, single reads for ONT'
                       ))
    parser.add_argument('--host_fraction', default="RANDOM",
                       help='Fraction of host reads')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility')
    parser.add_argument('--model', default="hiseq",
                       help='Sequencing model to use')
    parser.add_argument('--cores', type=int, default=1,
                       help='Number of cores for snakemake')
    parser.add_argument(
        '--simulator',
        default='iss',
        choices=['iss', 'camisim'],
        help='Read/community generator: InSilicoSeq (default) or optional CAMISIM',
    )
    parser.add_argument(
        '--reads_generator',
        '--reads-generator',
        dest='reads_generator',
        default=None,
        help=(
            'iss|art|wgsim|nanosim|hybrid|camisim, or a name imported with '
            '`samovar tools import --type reads`'
        ),
    )
    parser.add_argument(
        '-i',
        '--abundance',
        dest='abundance_table',
        default=None,
        help='Optional abundance CSV (taxid, N_<sample>…) instead of a uniform community',
    )
    parser.add_argument(
        '--flags',
        nargs=2,
        action='append',
        metavar=('TOOL', 'FLAGS'),
        dest='tool_flags',
        default=None,
        help='Extra flags for a reads_generator, e.g. --flags iss "--gc_bias"',
    )
    parser.add_argument(
        '--camisim-mode',
        default=None,
        help='CAMISIM mode: table (abundance then ISS), illumina, ont, wgsim, hybrid',
    )
    parser.add_argument(
        '--camisim-config',
        default=None,
        help='Existing CAMISIM YAML to overlay (still overwritten by CLI critical values)',
    )
    parser.add_argument(
        '--size-gbp',
        type=float,
        default=None,
        help='CAMISIM sample size in Gbp (default: derived from --total_reads)',
    )
    
    return parser.parse_args(argv)

def setup_iss_test(args: Optional[argparse.Namespace] = None) -> Dict[str, str]:
    """Main function to set up the ISS or CAMISIM generate configuration"""
    if args is None:
        args = parse_args()

    from samovar.main_config import merge_flag_strings
    from samovar.reads_generators import (
        attach_reads_flags,
        camisim_mode_for_reads_generator,
        flags_apply_to_reads_generator,
        require_known_reads_generator,
        resolve_reads_generator,
        write_custom_generate_script,
    )

    raw_name = getattr(args, "reads_generator", None)
    simulator = str(getattr(args, "simulator", None) or "iss").strip().lower()
    camisim_mode = getattr(args, "camisim_mode", None)
    if not raw_name:
        if camisim_mode and simulator == "iss":
            simulator = "camisim"
        if simulator == "camisim":
            mode = str(camisim_mode or "table").strip().lower()
            if mode in {"", "table", "community", "abundance", "design"}:
                raw_name = "camisim"
            elif mode in {"hybrid", "mix", "multi"}:
                raw_name = "hybrid"
            else:
                raw_name = mode
        else:
            raw_name = "iss"
    kind, canon = resolve_reads_generator(raw_name)
    if kind == "custom":
        require_known_reads_generator(canon)
    args.reads_generator = canon

    flag_parts = [getattr(args, "extra_flags", None)]
    for item in getattr(args, "tool_flags", None) or []:
        if not item or len(item) < 2:
            continue
        target, flags = item[0], item[1]
        if flags_apply_to_reads_generator(target, canon):
            flag_parts.append(flags)
    merged = attach_reads_flags(canon, {"extra_flags": merge_flag_strings(*flag_parts)})
    args.extra_flags = merged.get("extra_flags") or ""

    mapped_mode = camisim_mode_for_reads_generator(canon, camisim_mode)
    if mapped_mode is not None:
        from samovar.camisim import setup_camisim_generate

        args.camisim_mode = mapped_mode
        if getattr(args, "abundance_table", None):
            # CAMISIM generate still uses genome_dir; abundance is stored on the YAML.
            pass
        result = setup_camisim_generate(args)
        yaml_path = result.get("config") or result.get("yaml")
        if yaml_path and args.extra_flags:
            from pathlib import Path
            import yaml

            data = yaml.safe_load(Path(yaml_path).read_text(encoding="utf-8")) or {}
            data["extra_flags"] = args.extra_flags
            data["reads_generator_flags"] = args.extra_flags
            data["reads_generator"] = canon
            if getattr(args, "abundance_table", None):
                from samovar.paths import absolute_path

                data["abundance_table"] = absolute_path(args.abundance_table)
            Path(yaml_path).write_text(yaml.dump(data), encoding="utf-8")
        return result

    accessions = list(getattr(args, "accessions", None) or [])
    reindex_mode = int(getattr(args, "reindex", 0) or 0)
    keep_raw = bool(int(getattr(args, "raw_genomes", 0) or 0))
    if accessions:
        from samovar.genome_fetcher import default_entrez_email, materialize_accessions
        from samovar.genome_index import run_processed_dir

        paths = materialize_accessions(
            accessions,
            output_dir=args.output_dir,
            email=default_entrez_email(),
            reindex_mode=reindex_mode,
            keep_raw=keep_raw,
        )
        if not paths:
            raise SystemExit("generate: no genomes materialized from --accessions")
        args.genome_dir = str(run_processed_dir(args.output_dir))
    elif not getattr(args, "genome_dir", None):
        raise SystemExit("generate: --genome_dir or --accessions is required")

    config = ISSTestConfig.from_args(args)
    config_path = config.generate_config(config.output_dir)
    if kind == "custom":
        pipeline_path = write_custom_generate_script(
            config.output_dir, config_path, cores=config.cores
        )
    else:
        pipeline_path = config.generate_pipeline(config.output_dir)

    return {
        'config': config_path,
        'pipeline': pipeline_path
    } 