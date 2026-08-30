import os
import yaml
import argparse
from typing import Any, Dict, List, Optional, Sequence
from dataclasses import dataclass
from pathlib import Path
import sys

from samovar.genome_fetcher import argparse_max_genome_mb
from samovar.paths import add_output_dir_argument
from samovar.regenerate import UNLIMITED_MAX_GENOMES, argparse_max_genomes, parse_max_genomes


def parse_genome_tokens(values: Optional[Sequence[Any]]) -> List[str]:
    tokens: List[str] = []
    for raw in values or []:
        if raw in (None, False, ""):
            continue
        if isinstance(raw, (list, tuple)):
            tokens.extend(parse_genome_tokens(list(raw)))
            continue
        for piece in str(raw).replace(",", " ").split():
            if piece.strip():
                tokens.append(piece.strip())
    return tokens


def list_generate_genome_paths(genome_dir: str) -> List[Path]:
    """FASTAs under ``genome_dir/processed`` then ``genome_dir`` (stable sort)."""
    from samovar.seqio import is_fasta_name

    root = Path(genome_dir).expanduser() if genome_dir else Path()
    found: List[Path] = []
    seen = set()
    if not root.is_dir():
        return found
    for folder in (root / "processed", root):
        if not folder.is_dir():
            continue
        for path in sorted(folder.iterdir()):
            if not path.is_file():
                continue
            if not is_fasta_name(path.name, protein=False, nucleotide=True):
                continue
            try:
                key = str(path.resolve())
            except OSError:
                key = str(path)
            if key in seen:
                continue
            seen.add(key)
            found.append(path)
    return found


def resolve_locked_genome_paths(genome_dir: str, tokens: Sequence[str]) -> List[Path]:
    from samovar.seqio import taxid_from_fasta_name

    listed = list_generate_genome_paths(genome_dir)
    by_name = {p.name: p for p in listed}
    by_stem: Dict[str, Path] = {}
    for path in listed:
        taxid = taxid_from_fasta_name(path) or path.stem.split(".")[0]
        by_stem.setdefault(str(taxid), path)
        by_stem.setdefault(path.stem, path)
    out: List[Path] = []
    seen = set()
    for tok in tokens:
        path = Path(tok).expanduser()
        hit: Optional[Path] = None
        if path.is_file():
            hit = path
        elif tok in by_name:
            hit = by_name[tok]
        elif tok in by_stem:
            hit = by_stem[tok]
        if hit is None:
            continue
        try:
            key = str(hit.resolve())
        except OSError:
            key = str(hit)
        if key in seen:
            continue
        seen.add(key)
        out.append(hit)
    return out


def select_generate_genome_paths(
    genome_dir: str,
    max_genomes: Any = None,
    locked: Optional[Sequence[str]] = None,
) -> List[Path]:
    from samovar.regenerate import cap_sequence

    tokens = parse_genome_tokens(locked)
    if tokens:
        paths = resolve_locked_genome_paths(genome_dir, tokens)
        if paths:
            return paths
    return cap_sequence(list_generate_genome_paths(genome_dir), max_genomes)


def genome_lock_records(
    paths: Sequence[Any],
    host_genome: str = "",
) -> List[Dict[str, Any]]:
    from samovar.seqio import taxid_from_fasta_name

    host_key = ""
    if host_genome:
        try:
            host_key = str(Path(host_genome).expanduser().resolve())
        except OSError:
            host_key = str(Path(host_genome).expanduser())
    records: List[Dict[str, Any]] = []
    seen = set()
    for raw in paths:
        path = Path(raw)
        try:
            key = str(path.resolve()) if path.exists() else str(path)
        except OSError:
            key = str(path)
        if key in seen:
            continue
        seen.add(key)
        taxid = taxid_from_fasta_name(path) or path.stem.split(".")[0]
        records.append(
            {
                "taxid": str(taxid),
                "file": path.name,
                "path": key,
                "host": bool(host_key and key == host_key),
            }
        )
    if host_genome:
        host_path = Path(host_genome).expanduser()
        if host_path.is_file() and not any(rec.get("host") for rec in records):
            taxid = taxid_from_fasta_name(host_path) or host_path.stem.split(".")[0]
            records.append(
                {
                    "taxid": str(taxid),
                    "file": host_path.name,
                    "path": host_key or str(host_path),
                    "host": True,
                }
            )
    return records


def _genomes_from_generate_yaml(config_path: Optional[str]) -> Dict[str, Any]:
    if not config_path:
        return {}
    path = Path(config_path)
    if not path.is_file():
        return {}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    if not isinstance(data, dict):
        return {}
    files = data.get("genomes") or []
    used = genome_lock_records(files, str(data.get("host_genome") or ""))
    payload: Dict[str, Any] = {"used": used}
    if "max_genomes" in data:
        payload["max_genomes"] = data.get("max_genomes")
    return payload


def _finish_generate(args, result):
    try:
        from samovar.repro import _argv_from_args, _ns_to_dict, record_stage

        extra = {}
        locked = _genomes_from_generate_yaml((result or {}).get("config"))
        if locked.get("used"):
            extra["genomes"] = locked
        parsed = _ns_to_dict(args)
        argv = sys.argv[1:]
        generate_flags = ("--genome_dir", "--genome-dir", "--accessions", "--output_dir")
        if parsed and not any(
            str(tok).split("=", 1)[0] in generate_flags for tok in argv
        ):
            argv = _argv_from_args("generate", parsed)
        record_stage(
            "generate",
            getattr(args, "output_dir", None),
            args=args,
            argv=argv,
            extra=extra or None,
        )
    except Exception as exc:
        print(f"Warning: could not write Hydra snapshot: {exc}", file=sys.stderr)
    return result


def _materialize_generate_accessions(args):
    from samovar.genome_fetcher import default_entrez_email, materialize_accessions
    from samovar.genome_index import run_processed_dir

    accessions = list(getattr(args, "accessions", None) or [])
    if not accessions:
        return
    reindex_mode = int(getattr(args, "reindex", 0) or 0)
    keep_raw = bool(int(getattr(args, "raw_genomes", 0) or 0))
    paths = materialize_accessions(
        accessions,
        output_dir=args.output_dir,
        email=default_entrez_email(),
        reindex_mode=reindex_mode,
        keep_raw=keep_raw,
        max_genome_mb=getattr(args, "max_genome_mb", None),
        genome_skip_list=getattr(args, "genome_skip_list", None),
    )
    if not paths:
        raise SystemExit("generate: no genomes materialized from --accessions")
    args.genome_dir = str(run_processed_dir(args.output_dir))


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
    metagenome_generator: str = ""
    abundance_table: str = ""
    max_genome_mb: float = float("inf")
    genome_skip_list: str = ""
    max_genomes: float = UNLIMITED_MAX_GENOMES
    locked_genomes: Optional[List[str]] = None

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'ISSTestConfig':
        from samovar.paths import absolute_path

        abundance = getattr(args, "abundance_table", None) or ""
        if abundance:
            abundance = absolute_path(abundance)
        extra = getattr(args, "extra_flags", None) or ""
        from samovar.genome_fetcher import parse_genome_skip_list, parse_max_genome_mb

        skip = getattr(args, "genome_skip_list", None)
        skip_text = ""
        if skip is not None:
            skip_text = ",".join(sorted(parse_genome_skip_list(skip, default_from_env=False)))
        max_mb = parse_max_genome_mb(
            getattr(args, "max_genome_mb", None),
            default_from_env=getattr(args, "max_genome_mb", None) is None,
        )
        max_g = parse_max_genomes(
            getattr(args, "max_genomes", None),
            default_from_env=getattr(args, "max_genomes", None) is None,
        )
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
            metagenome_generator=str(getattr(args, "metagenome_generator", None) or ""),
            abundance_table=str(abundance or ""),
            max_genome_mb=max_mb,
            genome_skip_list=skip_text,
            max_genomes=max_g,
            locked_genomes=parse_genome_tokens(getattr(args, "genomes", None)),
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
            'genomes': [
                str(path)
                for path in select_generate_genome_paths(
                    self.genome_dir,
                    self.max_genomes,
                    self.locked_genomes,
                )
            ],
            'reads_generator': self.reads_generator,
        }
        if self.metagenome_generator:
            config['metagenome_generator'] = self.metagenome_generator
        if self.extra_flags:
            config['extra_flags'] = self.extra_flags
            config['reads_generator_flags'] = self.extra_flags
            if self.metagenome_generator:
                config['metagenome_generator_flags'] = self.extra_flags
        if self.abundance_table:
            config['abundance_table'] = self.abundance_table
        config['max_genome_mb'] = float(self.max_genome_mb)
        config['max_genomes'] = float(self.max_genomes)
        if self.genome_skip_list:
            config['genome_skip_list'] = self.genome_skip_list
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        return str(config_path)

    def generate_pipeline(self, base_dir: str) -> str:
        """Generate the ISS test pipeline script and return its path"""
        from samovar.paths import absolute_path, python_path, repo_root, runtime_path_prefix
        from samovar.paths import shell_outdir_override_snippet

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
{shell_outdir_override_snippet()}
mkdir -p "$out_dir"

# Generate reads with InSilicoSeq
snakemake -s {snakefile} \\
    --configfile {config_path} \\
    --directory "$out_dir" \\
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
    add_output_dir_argument(
        parser,
        required=True,
        help='Output directory for generated files',
    )
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
        '--threads',
        type=int,
        default=None,
        help='Alias of --cores; mapped to the generator via flags-translate (--cpus for ISS)',
    )
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
        '--metagenome_generator',
        '--metagenome-generator',
        dest='metagenome_generator',
        default=None,
        help=(
            'Combined community+reads tool: camisim|hybrid|nanosim, or a name '
            'imported with `samovar tools import --type meta`'
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
        help='Extra flags for a generator, e.g. --flags iss "--gc_bias" or --flags metagenome_generator "--model hiseq"',
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
    parser.add_argument(
        '--max-genomes', '--max_genomes', '--top',
        dest='max_genomes',
        type=argparse_max_genomes,
        default=None,
        help='Keep this many genomes from genome_dir (default: inf). The chosen files are written to Hydra.',
    )
    parser.add_argument(
        '--genomes',
        action='append',
        default=None,
        help='Lock this run to these FASTA paths/taxids (repeatable or comma-separated). Overrides --top scan order.',
    )
    parser.add_argument(
        '--max-genome-mb',
        dest='max_genome_mb',
        type=argparse_max_genome_mb,
        default=None,
        help='Skip new NCBI downloads larger than this many MB (default: inf). Cached genomes are kept.',
    )
    parser.add_argument(
        '--genome-skip-list',
        dest='genome_skip_list',
        default=None,
        help='Comma-separated taxids/accessions not to download (already-cached copies are still used)',
    )
    from samovar.tool_spec import apply_translated_flags, parse_known_with_custom

    args, leftover = parse_known_with_custom(parser, argv)
    if leftover:
        parser.error("unrecognized arguments: " + " ".join(leftover))
    if getattr(args, "threads", None):
        args.cores = int(args.threads)
    custom = dict(getattr(args, "custom_cli_flags", None) or {})
    if custom.get("--threads"):
        try:
            args.cores = int(custom["--threads"])
        except (TypeError, ValueError):
            pass
    gen = str(getattr(args, "reads_generator", None) or args.simulator or "iss")
    canon = {"--threads": args.cores, "--cores": args.cores, **custom}
    if str(gen).lower() in {"iss", "insilicoseq"}:
        canon.pop("--threads", None)
        canon.pop("--cores", None)
    extra = apply_translated_flags(
        str(getattr(args, "extra_flags", None) or ""),
        name=gen,
        canonical=canon,
    )
    args.extra_flags = extra
    return args

def setup_iss_test(args: Optional[argparse.Namespace] = None) -> Dict[str, str]:
    """Main function to set up the ISS or CAMISIM generate configuration"""
    if args is None:
        args = parse_args()

    from samovar.main_config import merge_flag_strings
    from samovar.metagenome_generators import (
        attach_metagenome_flags,
        flags_apply_to_metagenome_generator,
        require_known_metagenome_generator,
        resolve_metagenome_generator,
    )
    from samovar.reads_generators import (
        attach_reads_flags,
        camisim_mode_for_reads_generator,
        flags_apply_to_reads_generator,
        require_known_reads_generator,
        resolve_reads_generator,
        write_custom_generate_script,
    )

    meta_raw = getattr(args, "metagenome_generator", None)
    if meta_raw:
        mkind, mcanon = resolve_metagenome_generator(meta_raw)
        if mkind == "custom":
            require_known_metagenome_generator(mcanon)
        args.metagenome_generator = mcanon
        flag_parts = [getattr(args, "extra_flags", None)]
        for item in getattr(args, "tool_flags", None) or []:
            if not item or len(item) < 2:
                continue
            target, flags = item[0], item[1]
            if flags_apply_to_metagenome_generator(target, mcanon):
                flag_parts.append(flags)
        merged = attach_metagenome_flags(
            mcanon, {"extra_flags": merge_flag_strings(*flag_parts)}
        )
        args.extra_flags = merged.get("extra_flags") or ""
        mapped_mode = camisim_mode_for_reads_generator(mcanon, getattr(args, "camisim_mode", None))
        if mkind == "builtin" and mapped_mode is not None:
            from samovar.camisim import setup_camisim_generate

            args.camisim_mode = mapped_mode
            result = setup_camisim_generate(args)
            yaml_path = result.get("config") or result.get("yaml")
            if yaml_path:
                from pathlib import Path
                import yaml

                data = yaml.safe_load(Path(yaml_path).read_text(encoding="utf-8")) or {}
                if args.extra_flags:
                    data["extra_flags"] = args.extra_flags
                    data["metagenome_generator_flags"] = args.extra_flags
                data["metagenome_generator"] = mcanon
                Path(yaml_path).write_text(yaml.dump(data), encoding="utf-8")
            return _finish_generate(args, result)
        _materialize_generate_accessions(args)
        if not getattr(args, "genome_dir", None):
            raise SystemExit("generate: --genome_dir or --accessions is required")
        config = ISSTestConfig.from_args(args)
        config_path = config.generate_config(config.output_dir)
        pipeline_path = write_custom_generate_script(
            config.output_dir, config_path, cores=config.cores
        )
        return _finish_generate(args, {"config": config_path, "pipeline": pipeline_path})

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
        return _finish_generate(args, result)

    _materialize_generate_accessions(args)
    if not getattr(args, "genome_dir", None):
        raise SystemExit("generate: --genome_dir or --accessions is required")

    config = ISSTestConfig.from_args(args)
    config_path = config.generate_config(config.output_dir)
    if kind == "custom":
        pipeline_path = write_custom_generate_script(
            config.output_dir, config_path, cores=config.cores
        )
    else:
        pipeline_path = config.generate_pipeline(config.output_dir)

    return _finish_generate(
        args,
        {
            "config": config_path,
            "pipeline": pipeline_path,
        },
    ) 