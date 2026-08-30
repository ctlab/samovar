import os
import shutil
import sys
import yaml
import argparse
from typing import Dict, List, Optional, Union
from dataclasses import dataclass
from pathlib import Path

from samovar.exec_control import CHECKPOINT_STEPS, resolve_window
from samovar.paths import (
    absolute_path,
    add_output_dir_argument,
    ncbi_email,
    python_path,
    repo_root,
    runtime_path_prefix,
    shell_outdir_override_snippet,
)
from samovar.table_regenerators import flags_apply_to_regenerator, canonical_regeneration_modes
from samovar.table_scorers import flags_apply_to_table_scorer
from samovar.reads_generators import (
    flags_apply_to_reads_generator,
    require_known_reads_generator,
)
from samovar.metagenome_generators import (
    flags_apply_to_metagenome_generator,
    require_known_metagenome_generator,
)
from samovar.scorers import is_scoring_flag_target
from samovar.reprofilers import (
    REPROFILER_FLAG_GROUPS,
    is_reprofiler_flag_target,
    require_known_reprofiler,
)
from samovar.genome_resolve import normalize_reannotation_level
from samovar.genome_fetcher import (
    UNLIMITED_GENOME_MB,
    parse_genome_skip_list,
    parse_max_genome_mb,
)
from samovar.regenerate import UNLIMITED_MAX_GENOMES, parse_max_genomes
from samovar.ground_truth import (
    GROUND_TRUTH_TABLE,
    normalize_regenerated_mode,
)
from samovar.main_config import (
    flags_target_matches,
    imported_flags_for_names,
    iter_tools,
    merge_flag_strings,
)

DEFAULT_OUTPUT_DIR = "samovar_out"


def _default_email() -> str:
    return ncbi_email()


def _cmd_basename(cmd: str) -> str:
    if not cmd:
        return ""
    return Path(str(cmd).split()[0]).name.split(".")[0]


def _merge_annotator_launch_flags(
    annotators: List["AnnotatorConfig"],
    tool_flag_pairs: Optional[List],
) -> None:
    """Import ``tools.*[4]`` plus prepare ``--flags`` into each annotator's extra."""
    from samovar.paths import load_config

    tools = iter_tools(load_config())
    pairs = tool_flag_pairs or []
    groups = ("annotator", "annotators", "a", "ann")
    for ann in annotators:
        cmd_base = _cmd_basename(ann.cmd)
        imported = imported_flags_for_names(tools, ann.run_name, ann.type, cmd_base)
        launch = []
        for item in pairs:
            if not item or len(item) < 2:
                continue
            target, flags = item[0], item[1]
            if flags_target_matches(
                target, ann.run_name, ann.type, cmd_base, groups=groups
            ):
                launch.append(flags)
        merged = merge_flag_strings(imported, ann.extra, *launch)
        ann.extra = merged or None


def _apply_translated_cli_flags(config: "PipelineConfig", args: argparse.Namespace) -> None:
    """Push --threads/--cores (and custom-flags) through flags-translate into child extras."""
    from samovar.main_config import lookup_tool_record
    from samovar.paths import load_config
    from samovar.tool_spec import apply_translated_flags

    canonical = {
        "--threads": config.cores,
        "--cores": config.cores,
    }
    canonical.update(dict(getattr(args, "custom_cli_flags", None) or {}))
    cfg = load_config()
    native = {"kraken2", "kaiju", "kraken", "krakenuniq", "metaphlan", "metaphlan4"}
    for ann in config.annotators or []:
        name = ann.type or ann.run_name
        rec = lookup_tool_record(cfg, name) or lookup_tool_record(cfg, ann.run_name)
        extra_canon = dict(canonical)
        if str(ann.type or "").lower() in native:
            extra_canon.pop("--threads", None)
            extra_canon.pop("--cores", None)
        ann.extra = apply_translated_flags(
            ann.extra or "",
            name=name,
            record=rec,
            canonical=extra_canon,
        ) or None
    config.reads_generator_flags = apply_translated_flags(
        config.reads_generator_flags or "",
        name=config.reads_generator or "iss",
        record=lookup_tool_record(cfg, config.reads_generator or "iss"),
        canonical=(
            {k: v for k, v in canonical.items() if k not in {"--threads", "--cores"}}
            if str(config.reads_generator or "iss").lower() in {"iss", "insilicoseq"}
            else canonical
        ),
    ) or None
    if config.metagenome_generator:
        config.metagenome_generator_flags = apply_translated_flags(
            config.metagenome_generator_flags or "",
            name=config.metagenome_generator,
            record=lookup_tool_record(cfg, config.metagenome_generator),
            canonical=canonical,
        ) or None
    config.regeneration_extra_flags = apply_translated_flags(
        config.regeneration_extra_flags or "",
        name=(config.regeneration_mode or "direct"),
        record=lookup_tool_record(cfg, str(config.regeneration_mode or "")),
        canonical=canonical,
    ) or None
    config.scoring_flags = apply_translated_flags(
        config.scoring_flags or "",
        name="scoring",
        canonical=canonical,
    ) or None
    named_s = dict(config.scoring_tool_flags or {})
    for sname in list(named_s) + list(config.scoring_tools or []):
        named_s[sname] = apply_translated_flags(
            named_s.get(sname) or "",
            name=sname,
            record=lookup_tool_record(cfg, sname),
            canonical=canonical,
        )
    config.scoring_tool_flags = named_s or {}
    config.reprofiler_flags = apply_translated_flags(
        config.reprofiler_flags or "",
        name=config.reprofiler or "ensemble",
        record=lookup_tool_record(cfg, str(config.reprofiler or "")),
        canonical=canonical,
    ) or None


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


def _expand_indexed_database(cmd: str, db_path: str, extra: Optional[str]) -> tuple:
    """If ``db_path`` is a registered database name, substitute path and stored flags."""
    from samovar.db_spec import annotator_flags_for_db, lookup_database_record
    from samovar.paths import load_config

    token = (db_path or "").strip()
    if not token or token == ".":
        return db_path, extra
    as_path = Path(token).expanduser()
    if as_path.exists():
        return db_path, extra
    tool = Path(cmd).name.split(".")[0]
    rec = lookup_database_record(load_config(), tool, token)
    if rec is None:
        return db_path, extra
    resolved = str(rec.get("path") or "").strip() or db_path
    stored = annotator_flags_for_db(tool, resolved, str(rec.get("flags") or ""))
    merged = " ".join(p for p in (stored, extra or "") if p).strip() or None
    return resolved, merged


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
    max_genomes: float = UNLIMITED_MAX_GENOMES
    regeneration_mode: str = "direct"
    regeneration_modes: Optional[List[str]] = None
    table_score: str = "shannon_ks"
    regeneration_n: Optional[int] = None
    regeneration_n_reads: int = 1000
    regeneration_seed: int = 42
    reannotation_level: str = "taxid"
    rescale_abundance: bool = False
    samples_metadata: Optional[str] = None
    regeneration_extra_flags: Optional[str] = None
    reads_generator: str = "iss"
    reads_generator_flags: Optional[str] = None
    metagenome_generator: Optional[str] = None
    metagenome_generator_flags: Optional[str] = None
    gzip_genomes: bool = True
    gzip_reads: bool = False
    reuse_genomes: bool = True
    use_test_genomes: bool = False
    genome_dirs: Optional[List[str]] = None
    max_genome_mb: float = UNLIMITED_GENOME_MB
    genome_skip_list: Optional[List[str]] = None
    run_multiqc: Optional[bool] = None
    initial_ground_truth_table: Optional[str] = None
    regenerated_metagenomes: str = "parse-genome"
    scoring_flags: Optional[str] = None
    scoring_tools: Optional[List[str]] = None
    scoring_tool_flags: Optional[Dict[str, str]] = None
    reprofiler: str = "ensemble"
    reprofiler_flags: Optional[str] = None
    reprofiler_tool_flags: Optional[Dict[str, str]] = None
    startpoint: str = "setup_reads"
    endpoint: str = "viz_reprofiled"
    export_formats: Optional[List[str]] = None
    export_taxdump: Optional[str] = None
    export_taxonomy: Optional[str] = None

    def __post_init__(self):
        if self.annotators is None:
            self.annotators = []
        if self.genome_dirs is None:
            self.genome_dirs = []
        if self.genome_skip_list is None:
            self.genome_skip_list = []
        if self.scoring_tool_flags is None:
            self.scoring_tool_flags = {}
        if self.reprofiler_tool_flags is None:
            self.reprofiler_tool_flags = {}
        if not self.regeneration_modes:
            self.regeneration_modes = [self.regeneration_mode]
        if self.export_formats is None:
            self.export_formats = []
        if not self.email:
            self.email = _default_email()

    @classmethod
    def from_args(cls, args: argparse.Namespace) -> 'PipelineConfig':
        config = cls()
        config.cores = getattr(args, 'cores', 1) or 1
        threads = getattr(args, "threads", None)
        if threads:
            config.cores = int(threads)
        custom_cli = dict(getattr(args, "custom_cli_flags", None) or {})
        if custom_cli.get("--threads"):
            try:
                config.cores = int(custom_cli["--threads"])
            except (TypeError, ValueError):
                pass
        if custom_cli.get("--cores"):
            try:
                config.cores = int(custom_cli["--cores"])
            except (TypeError, ValueError):
                pass
        config.max_genomes = parse_max_genomes(
            getattr(args, "max_genomes", None),
            default_from_env=getattr(args, "max_genomes", None) is None,
        )
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
                yaml_mode = (
                    input_config.get("table_reads_generators")
                    or input_config.get("table_reads_generator")
                    or input_config.get("regeneration_mode")
                    or config.regeneration_mode
                )
                yaml_modes = canonical_regeneration_modes(yaml_mode) or ["direct"]
                config.regeneration_modes = yaml_modes
                config.regeneration_mode = yaml_modes[0]
                yaml_score = (
                    input_config.get("table_score")
                    or input_config.get("table_reads_scorer")
                    or input_config.get("table-score")
                )
                if yaml_score:
                    config.table_score = str(yaml_score)
                meta = input_config.get("samples_metadata") or input_config.get(
                    "metadata"
                )
                if meta not in (None, "", False):
                    config.samples_metadata = absolute_path(str(meta))
                yaml_flags = input_config.get("table_reads_generator_flags") or input_config.get(
                    "extra_flags"
                )
                if yaml_flags:
                    config.regeneration_extra_flags = str(yaml_flags)
                yaml_reads = (
                    input_config.get("reads_generator")
                    or input_config.get("reads-generator")
                    or config.reads_generator
                )
                config.reads_generator = require_known_reads_generator(yaml_reads)
                yaml_rflags = input_config.get("reads_generator_flags")
                if yaml_rflags:
                    config.reads_generator_flags = str(yaml_rflags)
                yaml_meta = input_config.get("metagenome_generator")
                if yaml_meta:
                    config.metagenome_generator = require_known_metagenome_generator(
                        yaml_meta
                    )
                yaml_mflags = input_config.get("metagenome_generator_flags")
                if yaml_mflags:
                    config.metagenome_generator_flags = str(yaml_mflags)
                yaml_sflags = input_config.get("scoring_flags")
                if yaml_sflags:
                    config.scoring_flags = str(yaml_sflags)
                yaml_stools = input_config.get("scoring_tools")
                if yaml_stools not in (None, "", False):
                    if isinstance(yaml_stools, (list, tuple)):
                        config.scoring_tools = [str(x) for x in yaml_stools]
                    else:
                        config.scoring_tools = [
                            p.strip()
                            for p in str(yaml_stools).replace(",", " ").split()
                            if p.strip()
                        ]
                yaml_stf = input_config.get("scoring_tool_flags")
                if isinstance(yaml_stf, dict):
                    config.scoring_tool_flags = {
                        str(k): str(v) for k, v in yaml_stf.items()
                    }
                yaml_reprofiler = (
                    input_config.get("reprofiler")
                    or input_config.get("ml")
                    or input_config.get("reprofiling")
                )
                if yaml_reprofiler:
                    config.reprofiler = require_known_reprofiler(yaml_reprofiler)
                yaml_rfl = input_config.get("reprofiler_flags") or input_config.get(
                    "ml_flags"
                )
                if yaml_rfl:
                    config.reprofiler_flags = str(yaml_rfl)
                yaml_rtf = input_config.get("reprofiler_tool_flags")
                if isinstance(yaml_rtf, dict):
                    config.reprofiler_tool_flags = {
                        str(k): str(v) for k, v in yaml_rtf.items()
                    }
                config.reannotation_level = normalize_reannotation_level(
                    input_config.get(
                        "reannotation_level",
                        input_config.get("reannotation-level", config.reannotation_level),
                    )
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
                yaml_max_g = input_config.get("max_genomes", input_config.get("max-genomes"))
                if yaml_max_g is not None:
                    config.max_genomes = parse_max_genomes(
                        yaml_max_g, default_from_env=False
                    )
                yaml_max_mb = input_config.get("max_genome_mb", input_config.get("max-genome-mb"))
                if yaml_max_mb is not None:
                    config.max_genome_mb = parse_max_genome_mb(
                        yaml_max_mb, default_from_env=False
                    )
                yaml_skip = input_config.get("genome_skip_list", input_config.get("genome-skip-list"))
                if yaml_skip is not None:
                    config.genome_skip_list = sorted(
                        parse_genome_skip_list(yaml_skip, default_from_env=False)
                    )
                yaml_start = (
                    input_config.get("startpoint")
                    or input_config.get("start_point")
                    or input_config.get("start")
                )
                if yaml_start:
                    config.startpoint = str(yaml_start)
                yaml_end = (
                    input_config.get("endpoint")
                    or input_config.get("end_point")
                    or input_config.get("end")
                )
                if yaml_end:
                    config.endpoint = str(yaml_end)
                yaml_gt = (
                    input_config.get("initial_ground_truth_table")
                    or input_config.get("initial-ground-truth-table")
                    or input_config.get("ground_truth_table")
                )
                if yaml_gt:
                    config.initial_ground_truth_table = str(yaml_gt)
                yaml_regen = (
                    input_config.get("regenerated_metagenomes")
                    or input_config.get("regenerated-metagenomes")
                )
                if yaml_regen:
                    config.regenerated_metagenomes = normalize_regenerated_mode(yaml_regen)
                yaml_export = (
                    input_config.get("export_formats")
                    or input_config.get("export-formats")
                    or input_config.get("to")
                )
                if yaml_export:
                    from samovar.annotation_convert import parse_export_formats

                    config.export_formats = parse_export_formats(yaml_export)
                yaml_taxdump = (
                    input_config.get("export_taxdump")
                    or input_config.get("taxdump")
                )
                if yaml_taxdump:
                    config.export_taxdump = str(yaml_taxdump)
                yaml_taxonomy = (
                    input_config.get("export_taxonomy")
                    or input_config.get("taxonomy")
                )
                if yaml_taxonomy:
                    from samovar.taxonomy import normalize_taxonomy

                    config.export_taxonomy = normalize_taxonomy(yaml_taxonomy)
                
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
                    # Second part is the database path or indexed name
                    db_path = parts[1] if len(parts) > 1 else "."
                    extra = ' '.join(parts[2:]) if len(parts) > 2 else None
                    db_path, extra = _expand_indexed_database(cmd, db_path, extra)
                    
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
                    # Second part is the database path or indexed name
                    db_path = parts[1] if len(parts) > 1 else "."
                    extra = ' '.join(parts[2:]) if len(parts) > 2 else None
                    db_path, extra = _expand_indexed_database(cmd, db_path, extra)
                    
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
        cli_gt = getattr(args, "initial_ground_truth_table", None)
        if cli_gt:
            config.initial_ground_truth_table = absolute_path(cli_gt)
        cli_regen = getattr(args, "regenerated_metagenomes", None)
        if cli_regen:
            config.regenerated_metagenomes = normalize_regenerated_mode(cli_regen)
        if config.initial_ground_truth_table:
            config.initial_ground_truth_table = absolute_path(
                config.initial_ground_truth_table
            )
        config.regenerated_metagenomes = normalize_regenerated_mode(
            config.regenerated_metagenomes
        )
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
        if getattr(args, "max_genomes", None) is not None:
            config.max_genomes = parse_max_genomes(
                args.max_genomes, default_from_env=False
            )
        if getattr(args, "max_genome_mb", None) is not None:
            config.max_genome_mb = parse_max_genome_mb(
                args.max_genome_mb, default_from_env=False
            )
        if getattr(args, "genome_skip_list", None) is not None:
            config.genome_skip_list = sorted(
                parse_genome_skip_list(args.genome_skip_list, default_from_env=False)
            )

        # CLI regeneration depth overrides YAML when both are present.
        n_reads = getattr(args, "N_reads", None)
        if n_reads is None:
            n_reads = getattr(args, "n_reads", None)
        if n_reads is not None:
            config.regeneration_n_reads = int(n_reads)
        regen_n = getattr(args, "regeneration_n", None)
        if regen_n is not None:
            config.regeneration_n = int(regen_n)

        level = getattr(args, "reannotation_level", None)
        if level:
            config.reannotation_level = normalize_reannotation_level(level)

        cli_mode = getattr(args, "table_reads_generator", None) or getattr(
            args, "regeneration_mode", None
        )
        if cli_mode:
            cli_modes = canonical_regeneration_modes(cli_mode)
            if cli_modes:
                config.regeneration_modes = cli_modes
                config.regeneration_mode = cli_modes[0]
        cli_score = getattr(args, "table_score", None)
        if cli_score:
            config.table_score = str(cli_score)
        cli_reads = getattr(args, "reads_generator", None)
        if cli_reads:
            config.reads_generator = require_known_reads_generator(cli_reads)
        cli_meta_gen = getattr(args, "metagenome_generator", None)
        if cli_meta_gen:
            config.metagenome_generator = require_known_metagenome_generator(
                cli_meta_gen
            )
        cli_scoring = getattr(args, "scoring", None)
        if cli_scoring:
            tokens = [str(x).strip() for x in cli_scoring if str(x).strip()]
            if any(t.lower() in {"none", "off", "false", "0"} for t in tokens):
                config.scoring_tools = []
            else:
                config.scoring_tools = tokens
        cli_reprofiler = (
            getattr(args, "reprofiler", None)
            or getattr(args, "ml", None)
        )
        if cli_reprofiler:
            config.reprofiler = require_known_reprofiler(cli_reprofiler)
        cli_meta = getattr(args, "samples_metadata", None)
        if cli_meta:
            config.samples_metadata = absolute_path(str(cli_meta))
        pairs = getattr(args, "tool_flags", None) or []
        _merge_annotator_launch_flags(config.annotators, pairs)
        flag_parts = [config.regeneration_extra_flags]
        regen_targets = list(config.regeneration_modes or [config.regeneration_mode])
        sd2_score = str(config.table_score or "").lower().replace("-", "_") in {
            "sparsedossa2_cv",
            "sparsedossa2",
            "sd2_cv",
            "fitcv",
        }
        for item in pairs:
            if not item or len(item) < 2:
                continue
            target, flags = item[0], item[1]
            if any(flags_apply_to_regenerator(target, mode) for mode in regen_targets):
                flag_parts.append(flags)
                continue
            if flags_apply_to_table_scorer(target, config.table_score):
                flag_parts.append(flags)
                continue
            low = str(target).lower().replace("-", "_")
            if "sparsedossa2" in low or low in {"sd2", "sd2_cv", "sd2_fit"}:
                if any(str(m).startswith("sparsedossa2") for m in regen_targets) or sd2_score:
                    flag_parts.append(flags)
        merged = merge_flag_strings(*flag_parts)
        config.regeneration_extra_flags = merged or None
        rflag_parts = [config.reads_generator_flags]
        for item in pairs:
            if not item or len(item) < 2:
                continue
            target, flags = item[0], item[1]
            if flags_apply_to_reads_generator(target, config.reads_generator):
                rflag_parts.append(flags)
        config.reads_generator_flags = merge_flag_strings(*rflag_parts) or None
        mflag_parts = [config.metagenome_generator_flags]
        for item in pairs:
            if not item or len(item) < 2:
                continue
            target, flags = item[0], item[1]
            if config.metagenome_generator and flags_apply_to_metagenome_generator(
                target, config.metagenome_generator
            ):
                mflag_parts.append(flags)
        config.metagenome_generator_flags = merge_flag_strings(*mflag_parts) or None
        sflag_parts = [config.scoring_flags]
        named_sflags = dict(config.scoring_tool_flags or {})
        for item in pairs:
            if not item or len(item) < 2:
                continue
            target, flags = item[0], item[1]
            if not is_scoring_flag_target(str(target)):
                continue
            if flags_target_matches(
                str(target),
                groups=(
                    "scoring",
                    "score",
                    "viz",
                    "visualization",
                    "visualisation",
                    "visualizations",
                    "visualisations",
                    "plots",
                ),
            ):
                sflag_parts.append(flags)
            else:
                named_sflags[str(target)] = merge_flag_strings(
                    named_sflags.get(str(target)), flags
                )
        config.scoring_flags = merge_flag_strings(*sflag_parts) or None
        config.scoring_tool_flags = named_sflags or {}
        rprof_parts = [config.reprofiler_flags]
        named_rprof = dict(config.reprofiler_tool_flags or {})
        for item in pairs:
            if not item or len(item) < 2:
                continue
            target, flags = item[0], item[1]
            if not is_reprofiler_flag_target(str(target)):
                continue
            if flags_target_matches(str(target), groups=REPROFILER_FLAG_GROUPS):
                rprof_parts.append(flags)
            else:
                named_rprof[str(target)] = merge_flag_strings(
                    named_rprof.get(str(target)), flags
                )
        config.reprofiler_flags = merge_flag_strings(*rprof_parts) or None
        config.reprofiler_tool_flags = named_rprof or {}

        cli_start = getattr(args, "startpoint", None)
        if cli_start:
            config.startpoint = str(cli_start)
        cli_end = getattr(args, "endpoint", None)
        if cli_end:
            config.endpoint = str(cli_end)
        config.startpoint, config.endpoint = resolve_window(
            config.startpoint, config.endpoint
        )
        from samovar.annotation_convert import formats_from_prepare_args, parse_export_formats

        cli_export = formats_from_prepare_args(args)
        if cli_export:
            config.export_formats = parse_export_formats(
                config.export_formats, cli_export
            )
        cli_taxdump = getattr(args, "export_taxdump", None) or getattr(args, "taxdump", None)
        if cli_taxdump:
            config.export_taxdump = str(cli_taxdump)
        cli_taxonomy = getattr(args, "export_taxonomy", None) or getattr(
            args, "taxonomy", None
        )
        if cli_taxonomy:
            from samovar.taxonomy import normalize_taxonomy

            config.export_taxonomy = normalize_taxonomy(cli_taxonomy)
        _apply_translated_cli_flags(config, args)
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
            'observed_abundance_dir': str(base_path / 'initial_abundance'),
            'abundance_dir': str(base_path / 'regenerated' / '.regenerated_abundance'),
            'genome_dir': str(base_path / 'genomes'),
            'output_dir': str(base_path / 'regenerated'),
            'email': self.email,
            'read_length': self.read_length,
            'coverage': self.coverage,
            'max_genomes': float(getattr(self, 'max_genomes', UNLIMITED_MAX_GENOMES)),
            'cores': self.cores,
            'regeneration_mode': self.regeneration_mode,
            'table_reads_generator': self.regeneration_mode,
            'table_score': self.table_score,
            'reannotation_level': self.reannotation_level,
            'N_reads': self.regeneration_n_reads,
            'seed': self.regeneration_seed,
            'rescale_abundance': self.rescale_abundance,
            'gzip_genomes': self.gzip_genomes,
            'gzip_reads': self.gzip_reads,
            'regenerated_metagenomes': normalize_regenerated_mode(self.regenerated_metagenomes),
            'max_genome_mb': float(getattr(self, "max_genome_mb", UNLIMITED_GENOME_MB)),
            'genome_skip_list': list(getattr(self, "genome_skip_list", None) or []),
            'reads_generator': self.reads_generator,
        }
        modes = list(self.regeneration_modes or [self.regeneration_mode])
        annotation2iss_config["table_reads_generators"] = modes
        if len(modes) > 1:
            annotation2iss_config["table_reads_generator"] = modes
        if self.metagenome_generator:
            annotation2iss_config['metagenome_generator'] = self.metagenome_generator
        if self.regeneration_n:
            annotation2iss_config['N'] = self.regeneration_n
        if self.samples_metadata:
            annotation2iss_config['samples_metadata'] = self.samples_metadata
        if self.regeneration_extra_flags:
            annotation2iss_config['table_reads_generator_flags'] = self.regeneration_extra_flags
        if self.reads_generator_flags:
            annotation2iss_config['reads_generator_flags'] = self.reads_generator_flags
        if self.metagenome_generator_flags:
            annotation2iss_config['metagenome_generator_flags'] = (
                self.metagenome_generator_flags
            )
        annotation2iss_path = configs_dir / 'config_annotation2iss.yaml'
        with open(annotation2iss_path, 'w') as f:
            yaml.dump(annotation2iss_config, f)
        configs['annotation2iss'] = str(annotation2iss_path)
        configs['abundance2iss'] = str(annotation2iss_path)

        scoring_config = {
            'output_dir': str(base_path),
            'scoring_flags': self.scoring_flags or "",
            'scoring_tool_flags': self.scoring_tool_flags or {},
        }
        if self.scoring_tools is not None:
            scoring_config['scoring_tools'] = self.scoring_tools
        scoring_path = configs_dir / 'config_scoring.yaml'
        with open(scoring_path, 'w') as f:
            yaml.dump(scoring_config, f)
        configs['scoring'] = str(scoring_path)

        reprofiling_config = {
            'output_dir': str(base_path / 'reprofiled_annotations'),
            'initial_dir': str(base_path / 'initial_annotations'),
            'regenerated_path': str(
                base_path / 'regenerated_annotations' / 'combined_annotation_table.csv'
            ),
            'ground_truth_dir': str(base_path / 'regenerated' / '.regenerated_abundance'),
            'reprofiler': self.reprofiler,
            'reprofiler_flags': self.reprofiler_flags or "",
            'reprofiler_tool_flags': self.reprofiler_tool_flags or {},
            'seed': self.regeneration_seed,
        }
        reprofiling_path = configs_dir / 'config_reprofiling.yaml'
        with open(reprofiling_path, 'w') as f:
            yaml.dump(reprofiling_config, f)
        configs['reprofiling'] = str(reprofiling_path)

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

    def _annotation_export_bash(self, src_rel: str, stage: str) -> str:
        """Convert combined annotation tables after a combine/reprofile step."""
        formats = list(self.export_formats or [])
        if not formats:
            return ""
        import shlex

        tos = " ".join("--to " + shlex.quote(fmt) for fmt in formats)
        tax = ""
        if self.export_taxdump:
            tax = " --taxdump " + shlex.quote(str(self.export_taxdump))
        if self.export_taxonomy:
            tax += " --taxonomy " + shlex.quote(str(self.export_taxonomy))
        return (
            f"\n$PYTHON_PATH -m samovar.annotation_convert \\\n"
            f'    -i "$out_dir/{src_rel}" -o "$out_dir/exports/{stage}"{tax} \\\n'
            f"    {tos}"
        )

    def generate_pipeline(self, base_dir: str) -> str:
        """Generate the pipeline script and return its path"""
        base_dir = absolute_path(base_dir)
        base_path = Path(base_dir)
        log_dir = base_path / '.log'
        log_dir.mkdir(parents=True, exist_ok=True)
        
        pipeline_path = log_dir / 'samovar.sh'
        start_s, end_s = resolve_window(self.startpoint, self.endpoint)
        (log_dir / "window.env").write_text(
            f"export SAMOVAR_START={start_s}\nexport SAMOVAR_END={end_s}\n",
            encoding="utf-8",
        )
        
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
        input_dir_export = self.input_dir or ""
        regen_mode = normalize_regenerated_mode(self.regenerated_metagenomes)
        inject_taxid = "0" if regen_mode == GROUND_TRUTH_TABLE else "1"
        gt_initial = ""
        if self.initial_ground_truth_table:
            gt_path = absolute_path(self.initial_ground_truth_table)
            gt_initial = f' \\\n    --truth-table "{gt_path}"'
        regen_truth = Path(base_dir) / ".log" / "regenerated_ground_truth.tsv"
        if regen_mode == GROUND_TRUTH_TABLE:
            combine_regen_body = f"""$PYTHON_PATH -m samovar.ground_truth from-genomes \\
    --genome-dir "$out_dir/genomes" --output "{regen_truth}" --overwrite
$PYTHON_PATH {wf / 'combine_annotation_tables.py'} \\
    -i "$out_dir/regenerated_reports" \\
    -o "$out_dir/regenerated_annotations" \\
    -s 2 \\
    --truth-table "{regen_truth}" """
        else:
            combine_regen_body = f"""$PYTHON_PATH {wf / 'combine_annotation_tables.py'} \\
    -i "$out_dir/regenerated_reports" \\
    -o "$out_dir/regenerated_annotations" \\
    -s 2"""
        combine_regen_body = combine_regen_body + self._annotation_export_bash(
            "regenerated_annotations", "regenerated"
        )
        
        # Generate pipeline script (absolute paths so exec works from any cwd).
        # Completed steps write $out_dir/.log/checkpoints/<name>.done and are
        # skipped on the next exec unless SAMOVAR_REDO=1 / --redo.
        # SAMOVAR_START / SAMOVAR_END select a slice; exec may override them.
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
{shell_outdir_override_snippet()}
CKPT="$out_dir/.log/checkpoints"
mkdir -p "$CKPT"
# Checkpoint names: {step_names}
export SAMOVAR_START="${{SAMOVAR_START:-{start_s}}}"
export SAMOVAR_END="${{SAMOVAR_END:-{end_s}}}"
export SAMOVAR_INPUT_DIR="${{SAMOVAR_INPUT_DIR:-{input_dir_export}}}"

ckpt_skip() {{
  if ! "$PYTHON_PATH" -m samovar.exec_control active "$1"; then
    return 0
  fi
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
mkdir -p "$out_dir/initial" "$out_dir/initial_reports" "$out_dir/initial_abundance" "$out_dir/regenerated" "$out_dir/regenerated_reports"
"$PYTHON_PATH" -m samovar.exec_control require "$out_dir"
# NCBI genome cache reuse (truncated data/test_genomes is not a library).
export SAMOVAR_REUSE_GENOMES="${{SAMOVAR_REUSE_GENOMES:-{reuse_flag}}}"
export SAMOVAR_ALLOW_TEST_GENOMES="${{SAMOVAR_ALLOW_TEST_GENOMES:-{test_flag}}}"
export SAMOVAR_MULTIQC="${{SAMOVAR_MULTIQC:-{multiqc_flag}}}"
export SAMOVAR_GENOME_DIRS="${{SAMOVAR_GENOME_DIRS:+$SAMOVAR_GENOME_DIRS:}}{extra_genome_dirs}"
export SAMOVAR_RUN_DIR="$out_dir"
export SAMOVAR_REGENERATED_METAGENOMES="{regen_mode}"
export SAMOVAR_INJECT_TAXID="${{SAMOVAR_INJECT_TAXID:-{inject_taxid}}}"
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
    --directory "$out_dir" \\
    --cores {self.cores}""",
        )

        combine_initial = _checkpoint_block(
            "combine_initial",
            f"""$PYTHON_PATH {wf / 'combine_annotation_tables.py'} \\
    -i "$out_dir/initial_reports" \\
    -o "$out_dir/initial_annotations"{gt_initial}"""
            + self._annotation_export_bash("initial_annotations", "initial"),
        )

        viz_initial = _checkpoint_block(
            "viz_initial",
            f"""$PYTHON_PATH {wf / 'compare_annotations.py'} \\
    --annotation_dir "$out_dir/initial_annotations" \\
    --output_dir "$out_dir/initial_annotations_plots" \\
    --show_top 0
"$PYTHON_PATH" -m samovar.scorers run --output_dir "$out_dir" --stage viz_initial \\
    --config {configs['scoring']}""",
        )

        abundance_tables = _checkpoint_block(
            "abundance_tables",
            """$PYTHON_PATH -m samovar.abundance materialize --output_dir "$out_dir\"""",
        )

        regenerate_tables = _checkpoint_block(
            "regenerate_tables",
            f"""$PYTHON_PATH -m samovar.regenerate tables \\
    --output_dir "$out_dir" \\
    --config {configs['annotation2iss']}""",
        )

        score_regenerated_tables = _checkpoint_block(
            "score_regenerated_tables",
            f"""$PYTHON_PATH -m samovar.table_scorers stage \\
    --output_dir "$out_dir" \\
    --config {configs['annotation2iss']}""",
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
            f"""snakemake -s {wf / 'abundance2iss' / 'Snakefile'} \\
    --configfile {configs['abundance2iss']} \\
    --directory "$out_dir" \\
    --cores {self.cores}
"$PYTHON_PATH" -m samovar.exec_control cleanup "$out_dir" || true""",
        )

        early_exit = """if "$PYTHON_PATH" -m samovar.exec_control needs-regen; then
if ! $PYTHON_PATH -c "from samovar.seqio import has_r1_reads; raise SystemExit(0 if has_r1_reads('$out_dir/regenerated') else 1)"; then
    echo "No regenerated reads were produced; skipping re-annotation and reprofiling."
    cleanup_tmp_if_requested
    exit 0
fi
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
    --directory "$out_dir" \\
    --cores {self.cores}""",
        )

        combine_regenerated = _checkpoint_block(
            "combine_regenerated",
            combine_regen_body,
        )

        viz_regenerated = _checkpoint_block(
            "viz_regenerated",
            f"""$PYTHON_PATH {wf / 'compare_annotations.py'} \\
    --annotation_dir "$out_dir/regenerated_annotations" \\
    --output_dir "$out_dir/regenerated_annotations_plots" \\
    --csv "$out_dir/regenerated_annotations/combined_annotation_table.csv" \\
    --show_top 0
"$PYTHON_PATH" -m samovar.scorers run --output_dir "$out_dir" --stage viz_regenerated \\
    --config {configs['scoring']}""",
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
    --ground-truth "$out_dir/regenerated/.regenerated_abundance" \\
    --output_dir "$out_dir/reprofiled_annotations" \\
    --reprofiler {self.reprofiler} \\
    --config {configs['reprofiling']} \\
    --seed {self.regeneration_seed} \\
    $FEATURE_ARG"""
            + self._annotation_export_bash("reprofiled_annotations", "reprofiled"),
        )

        viz_reprofiled = _checkpoint_block(
            "viz_reprofiled",
            f"""$PYTHON_PATH {wf / 'compare_annotations.py'} \\
    --annotation_dir "$out_dir/reprofiled_annotations" \\
    --output_dir "$out_dir/reprofiled_annotations_plots" \\
    --csv "$out_dir/reprofiled_annotations/combined_annotation_table.csv" \\
    --show_top 0
"$PYTHON_PATH" -m samovar.scorers run --output_dir "$out_dir" --stage viz_reprofiled \\
    --config {configs['scoring']}""",
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
                abundance_tables,
                regenerate_tables,
                score_regenerated_tables,
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
    add_output_dir_argument(
        parser,
        default=None,
        help='Output directory (default: samovar_out, or from --input_config)',
    )
    
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
    
    try:
        config = PipelineConfig.from_args(args)
    except ValueError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)
    try:
        from samovar.genome_cache import remember_prepare_genome_paths

        remember_prepare_genome_paths(extra_dirs=config.genome_dirs)
    except OSError as exc:
        print(f"Warning: could not update genome cache config: {exc}", file=sys.stderr)
    configs = config.generate_configs(config.output_dir)
    pipeline_path = config.generate_pipeline(config.output_dir)
    try:
        from samovar.repro import record_stage

        record_stage(
            "prepare",
            config.output_dir,
            args=args,
            argv=sys.argv[1:],
        )
    except Exception as exc:
        print(f"Warning: could not write Hydra snapshot: {exc}", file=sys.stderr)

    return {
        'configs': configs,
        'pipeline': pipeline_path
    } 