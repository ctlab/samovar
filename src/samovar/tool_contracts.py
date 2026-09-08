"""In→out contracts for ``samovar tools import`` groups.

Each group has a documented input and output. Pytest in
``tests/test_tool_contracts.py`` exercises a dest file against one group.
``samovar tools import --pytest`` runs that check before writing config.
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

from samovar.main_config import normalize_tool_group
from samovar.paths import repo_root


def _contract_repo_root() -> Path:
    """Repo that contains tests/test_tool_contracts.py (not the install config root)."""
    here = Path(__file__).resolve()
    candidate = here.parent.parent.parent
    if (candidate / "tests" / "test_tool_contracts.py").is_file():
        return candidate
    return repo_root()

CONTRACTS: Dict[str, Dict[str, str]] = {
    "annotator": {
        "in": "FASTQ R1/R2 via CLI -i/-I, -d DB, -o out, -t threads "
        "(or Python parse_output(path) on a two-column TSV).",
        "out": "TSV/table with columns seq, taxID (CustomAnnotator parse).",
    },
    "table_reads_generator": {
        "in": "Abundance or long annotation (Annotation / DataFrame / dir of "
        "taxid+N_* CSVs). Optional metadata DataFrame. config may include "
        "max_genomes (default inf).",
        "out": "dict[str, DataFrame] of taxid + N_<sample> tables from regenerate().",
    },
    "table_scoring": {
        "in": "Observed annotation/abundance plus generated tables "
        "(score_annotator(...) or score_table(observed, generated, config)).",
        "out": "dict with rank_value (and usually ok, pvalue, scorer).",
    },
    "scoring": {
        "in": "List of paths under the run dir (glob *annotations by default), "
        "output_dir, config.",
        "out": "Files written under dest from score(inputs, dest, config).",
    },
    "reads_generator": {
        "in": "job spec (abundance_table, genome_dir, output_dir, n_samples, …) "
        "plus optional metadata; generate(spec, metadata, config).",
        "out": "list[str] FASTQ paths (R1/R2).",
    },
    "metagenome_generator": {
        "in": "Same generate(spec, metadata, config) as reads_generator "
        "(community + sequencing in one tool). config may include max_genomes "
        "(default inf) to cap non-host genomes.",
        "out": "list[str] FASTQ paths.",
    },
    "reprofiler": {
        "in": "reprofile(regenerated_df, ground_truth_tables, initial_tables, config).",
        "out": "ReprofileResult or dict with tables keyed by sample.",
    },
    "annotation_converter": {
        "in": "Annotation (or Path when --from is custom) plus dest path and config "
        "(from/to, extra_argv). dump(annotation, dest, config) and/or load(path, config).",
        "out": "Files at dest (dump) or an Annotation (load). convert(src, dest, config) is also accepted.",
    },
    "export": {
        "in": "Annotation (per-read taxID_* / seq) plus dest and config. "
        "config may include reference (Annotation with taxID_true / true), "
        "to (abundance|kraken2|cami|…), extra_argv. "
        "export(annotation, dest, config).",
        "out": "Abundance-like tables at dest (taxid + N_<sample>, or kraken2/cami reports).",
    },
    "qc": {
        "in": "FASTQ R1 (and optional R2) plus dest paths and config "
        "(min_gc/max_gc, extra_argv). Python trim(r1, r2, dest_r1, dest_r2, config) "
        "or a native binary (fastp -i/-I/-o/-O, cutadapt -o/-p, "
        "trimmomatic PE/SE, chopper/nanofilt stdin→stdout).",
        "out": "Trimmed FASTQ at dest_r1/dest_r2 (same layout as input; empty files allowed).",
    },
    "assembler": {
        "in": "FASTQ R1/R2 via CLI -i/-I, -o contig FASTA, -t threads.",
        "out": "Contig FASTA at -o (empty file allowed if assembly produced none).",
    },
    "gene_caller": {
        "in": "Contig or MAG FASTA -c (file or directory), gene directory -o, -t threads.",
        "out": "Directory of amino-acid FASTA (*.faa) and GFF (*.gff); Prodigal default.",
    },
    "binner": {
        "in": "Contigs -c, optional reads -i/-I, MAG directory -o, -t threads.",
        "out": "Directory of MAG FASTA files (*.fa / *.fna).",
    },
    "binner_qc": {
        "in": "MAG directory -c, optional DB -d, QC table -o, -t threads.",
        "out": "TSV mag_id,completeness,contamination,score (plus native report).",
    },
    "binner_combine": {
        "in": "One or more MAG dirs -c dir1,dir2, optional QC table -q, -o MAG dir, -t.",
        "out": "One non-redundant MAG FASTA directory (best bins).",
    },
    "mag_taxonomy": {
        "in": "MAG directory -c, taxonomy DB -d, table -o, -t threads.",
        "out": "TSV mag_id,taxid,lineage.",
    },
    "aligner": {
        "in": "FASTQ -i/-I, MAG FASTA or dir -r, BAM -o, -t threads.",
        "out": "Sorted BAM at -o (index .bai when samtools is available).",
    },
    "mag_quantifier": {
        "in": "BAM -b, MAG dir -r, table -o.",
        "out": "MAG abundance table (mag_id plus count/coverage columns).",
    },
    "taxon_quantifier": {
        "in": "MAG abundance -a, MAG taxonomy -x, table -o.",
        "out": "taxid + N_* abundance table.",
    },
    "read_assigner": {
        "in": "BAM -b, MAG taxonomy -x, TSV -o.",
        "out": "Two-column seq, taxID TSV (Annotation-ready).",
    },
}

GROUP_TO_TESTNODE = {
    "annotator": "tests/test_tool_contracts.py::test_annotator_contract",
    "table_reads_generator": "tests/test_tool_contracts.py::test_table_regenerator_contract",
    "table_scoring": "tests/test_tool_contracts.py::test_table_scoring_contract",
    "scoring": "tests/test_tool_contracts.py::test_scoring_contract",
    "reads_generator": "tests/test_tool_contracts.py::test_reads_generator_contract",
    "metagenome_generator": "tests/test_tool_contracts.py::test_metagenome_generator_contract",
    "reprofiler": "tests/test_tool_contracts.py::test_reprofiler_contract",
    "annotation_converter": "tests/test_tool_contracts.py::test_annotation_converter_contract",
    "export": "tests/test_tool_contracts.py::test_export_contract",
    "qc": "tests/test_tool_contracts.py::test_qc_contract",
    "assembler": "tests/test_tool_contracts.py::test_assembler_contract",
    "gene_caller": "tests/test_tool_contracts.py::test_gene_caller_contract",
    "binner": "tests/test_tool_contracts.py::test_binner_contract",
    "binner_qc": "tests/test_tool_contracts.py::test_binner_qc_contract",
    "binner_combine": "tests/test_tool_contracts.py::test_binner_combine_contract",
    "mag_taxonomy": "tests/test_tool_contracts.py::test_mag_taxonomy_contract",
    "aligner": "tests/test_tool_contracts.py::test_aligner_contract",
    "mag_quantifier": "tests/test_tool_contracts.py::test_mag_quantifier_contract",
    "taxon_quantifier": "tests/test_tool_contracts.py::test_taxon_quantifier_contract",
    "read_assigner": "tests/test_tool_contracts.py::test_read_assigner_contract",
}

DEFAULT_TOOLS = {
    "annotator": "tests/tools/dummy_annotator.py",
    "table_reads_generator": "tests/tools/identity_table.py",
    "table_scoring": "tests/data/bray_ks_table_scorer.py",
    "scoring": "tests/tools/count_annotations.py",
    "reads_generator": "tests/tools/echo_reads.py",
    "metagenome_generator": "tests/tools/echo_reads.py",
    "reprofiler": "tests/tools/linear_wrapper.py",
    "annotation_converter": "tests/tools/echo_annotation_converter.py",
    "export": "tests/tools/identity_export.py",
    "qc": "tests/tools/gc_filter.py",
    "assembler": "tests/tools/dummy_assembler.py",
    "gene_caller": "tests/tools/dummy_gene_caller.py",
    "binner": "tests/tools/dummy_binner.py",
    "binner_qc": "tests/tools/dummy_binner_qc.py",
    "binner_combine": "tests/tools/dummy_binner_combine.py",
    "mag_taxonomy": "tests/tools/dummy_mag_taxonomy.py",
    "aligner": "tests/tools/dummy_aligner.py",
    "mag_quantifier": "tests/tools/dummy_mag_quantifier.py",
    "taxon_quantifier": "tests/tools/dummy_taxon_quantifier.py",
    "read_assigner": "tests/tools/dummy_read_assigner.py",
}


def format_contract(group: str) -> str:
    spec = CONTRACTS.get(group) or {}
    inn = spec.get("in", "(see wiki/pipeline.md)")
    out = spec.get("out", "")
    return f"type {group}\n  in:  {inn}\n  out: {out}"


def load_python_module(path: Path, name: str = "samovar_contract_tool"):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_contract_pytest(
    tool_path: str,
    tool_type: str,
    *,
    extra_args: Optional[list] = None,
) -> Tuple[int, str]:
    """Run the matching contract test. Return (returncode, combined output)."""
    group = normalize_tool_group(tool_type)
    if group not in GROUP_TO_TESTNODE:
        known = ", ".join(sorted(GROUP_TO_TESTNODE))
        raise ValueError(
            f"No in→out contract pytest for --type {tool_type} ({group}). "
            f"Known: {known}"
        )
    root = _contract_repo_root()
    node = GROUP_TO_TESTNODE[group]
    env = os.environ.copy()
    src = str(root / "src")
    env["PYTHONPATH"] = src + (":" + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    env["SAMOVAR_CONTRACT_TOOL"] = str(Path(tool_path).resolve())
    env["SAMOVAR_CONTRACT_TYPE"] = group
    cmd = [
        sys.executable,
        "-m",
        "pytest",
        str(root / node.split("::")[0]),
        "-q",
        "--tb=short",
        "-k",
        node.split("::")[1],
        "--tool",
        str(Path(tool_path).resolve()),
        "--tool-type",
        group,
    ]
    if extra_args:
        cmd.extend(extra_args)
    proc = subprocess.run(
        cmd,
        cwd=str(root),
        capture_output=True,
        text=True,
        env=env,
    )
    text = (proc.stdout or "") + (proc.stderr or "")
    return proc.returncode, text
