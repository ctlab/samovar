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

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'ISSTestConfig':
        return cls(
            genome_dir=args.genome_dir,
            output_dir=args.output_dir,
            host_genome=args.host_genome,
            n_samples=args.n_samples if args.n_samples is not None else 10,
            total_reads=args.total_reads if args.total_reads is not None else 2000,
            host_fraction=args.host_fraction if args.host_fraction is not None else "RANDOM",
            seed=args.seed if args.seed is not None else 42,
            model=args.model if args.model is not None else "hiseq"
        )

    def generate_config(self, base_dir: str) -> str:
        """Generate ISS test config file and return its path"""
        base_path = Path(base_dir)
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
            'genomes': []  # Will be automatically populated
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        return str(config_path)

    def generate_pipeline(self, base_dir: str) -> str:
        """Generate the ISS test pipeline script and return its path"""
        base_path = Path(base_dir)
        generate_dir = base_path / '.generate'
        generate_dir.mkdir(parents=True, exist_ok=True)
        
        pipeline_path = generate_dir / 'generate.sh'
        
        # Get config path
        config_path = self.generate_config(base_dir)
        
        # Generate pipeline script
        pipeline_content = f"""# Setup
set -e

if [ -f build/config.json ]; then
    PYTHON_PATH=$(grep -o '"python_path": *"[^"]*"' build/config.json | sed 's/"python_path": *"\\(.*\\)"/\\1/')
else
    echo "SamovaR is not installed: check build/config.json"
    exit 1
fi

out_dir="{base_dir}"
mkdir -p $out_dir

# Generate reads with InSilicoSeq
snakemake -s workflow/iss_test/Snakefile \\
    --configfile {config_path} \\
    --cores 1
"""
        
        with open(pipeline_path, 'w') as f:
            f.write(pipeline_content)
        
        return str(pipeline_path)

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='SamovaR ISS Test Configuration')
    
    # Required arguments
    parser.add_argument('--genome_dir', required=True,
                       help='Directory containing reference genomes')
    parser.add_argument('--output_dir', required=True,
                       help='Output directory for generated files')
    parser.add_argument('--host_genome', required=True,
                       help='Path to host genome file')
    
    # Optional arguments with defaults
    parser.add_argument('--n_samples', type=int, default=10,
                       help='Number of samples to generate')
    parser.add_argument('--total_reads', type=int, default=2000,
                       help='Total number of reads to generate')
    parser.add_argument('--host_fraction', default="RANDOM",
                       help='Fraction of host reads')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility')
    parser.add_argument('--model', default="hiseq",
                       help='Sequencing model to use')
    
    return parser.parse_args()

def setup_iss_test(args: Optional[argparse.Namespace] = None) -> Dict[str, str]:
    """Main function to set up the ISS test configuration"""
    if args is None:
        args = parse_args()
    
    config = ISSTestConfig.from_args(args)
    config_path = config.generate_config(config.output_dir)
    pipeline_path = config.generate_pipeline(config.output_dir)
    
    return {
        'config': config_path,
        'pipeline': pipeline_path
    } 