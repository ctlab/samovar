import glob
import os
import re
import shlex
import sqlite3
import sys
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

from samovar.baselines import CONSTANT_TAXID_NAMES
from samovar.paths import command_argv, python_path

DUMMY_TAXID = "9606"
DUMMY_TOOL_NAMES = CONSTANT_TAXID_NAMES
CONSTANT_TAXID_TOOL_NAMES = CONSTANT_TAXID_NAMES

# Thin CLI wrappers that already append ``-p <tool>`` and forward custom.sh flags.
CUSTOM_WRAPPER_SCRIPTS = {
    "custom.sh",
    "centrifuge.sh",
    "metax.sh",
    "metauto.sh",
    "assembly_hybrid.sh",
}

# These names always go through custom.sh (native CLIs are not the SamovaR flags).
CUSTOM_SH_ROUTER_TOOLS = {
    "centrifuge",
    "metax",
    "metauto",
    "assembly_hybrid",
} | DUMMY_TOOL_NAMES


def skip_empty_reads_cmd(input_r1: str, outputs: Sequence[str], run_cmd: str) -> str:
    """If R1 has no FASTQ records, touch expected outputs instead of running.

    Kraken2/Kaiju on a 0-byte or 1-byte FASTQ can exit 0 without writing ``--output``.
    """
    r1 = shlex.quote(str(input_r1))
    outs = " ".join(shlex.quote(str(path)) for path in outputs)
    py = shlex.quote(python_path())
    src_root = str(Path(__file__).resolve().parent.parent)
    code = (
        "from samovar.seqio import fastq_has_reads; import sys; "
        "raise SystemExit(0 if fastq_has_reads(sys.argv[1]) else 1)"
    )
    probe = (
        f"PYTHONPATH={shlex.quote(src_root)}${{PYTHONPATH:+:$PYTHONPATH}} "
        f"{py} -c {shlex.quote(code)} {r1}"
    )
    return f"if ! {probe}; then touch {outs}; else {run_cmd.strip()}; fi"


def _has_r2(input_r2: str) -> bool:
    return bool(str(input_r2 or "").strip())


def _r2_flag(flag: str, input_r2: str) -> str:
    """``-I path`` when paired; nothing when the Snakefile marked the sample single-end."""
    if not _has_r2(input_r2):
        return ""
    return f"{flag} {shlex.quote(str(input_r2))} "


def single_or_paired_reads_cmd(
    input_r1: str,
    input_r2: str,
    outputs: Sequence[str],
    *,
    single_cmd: str,
    paired_cmd: str,
) -> str:
    """Run a classifier in paired mode only when R2 contains reads.

    An empty ``input_r2`` is single-end: the paired command is not emitted.
    """
    if not _has_r2(input_r2):
        return skip_empty_reads_cmd(input_r1, outputs, single_cmd)
    r1 = shlex.quote(str(input_r1))
    r2 = shlex.quote(str(input_r2))
    outs = " ".join(shlex.quote(str(path)) for path in outputs)
    py = shlex.quote(python_path())
    src_root = str(Path(__file__).resolve().parent.parent)
    code = (
        "from samovar.seqio import fastq_has_reads; import sys; "
        "raise SystemExit(0 if fastq_has_reads(sys.argv[1]) else 1)"
    )
    prefix = f"PYTHONPATH={shlex.quote(src_root)}${{PYTHONPATH:+:$PYTHONPATH}} {py} -c {shlex.quote(code)}"
    probe_r1 = f"{prefix} {r1}"
    probe_r2 = f"{prefix} {r2}"
    return (
        f"if ! {probe_r1}; then touch {outs}; "
        f"elif {probe_r2}; then {paired_cmd.strip()}; "
        f"else {single_cmd.strip()}; fi"
    )


def _empty_taxid_frame() -> pd.DataFrame:
    return pd.DataFrame(columns=["seq", "taxID"])


def _read_table_or_empty(file_path: str) -> Optional[pd.DataFrame]:
    path = Path(file_path)
    try:
        if not path.is_file() or path.stat().st_size == 0:
            return None
    except OSError:
        return None
    try:
        df = pd.read_table(file_path, header=None)
    except (pd.errors.EmptyDataError, pd.errors.ParserError):
        return None
    if df is None or df.empty:
        return None
    return df


class BaseAnnotator:
    """Base class for taxonomic annotators.

    This class provides a unified interface for generating execution commands
    for Snakemake and parsing the resulting output files into a standardized format.
    """

    def __init__(self, run_config: Dict, config: Dict):
        """Initialize the annotator with run-specific and global configurations.

        Args:
            run_config: Dictionary containing settings for this specific run
                        (e.g., db_path, threads, extra flags).
            config: The global Snakemake configuration dictionary.
        """
        self.run_config = run_config
        self.config = config

        # Extract common parameters shared across all annotation tools
        self.run_name = run_config.get("run_name", "unknown_run")
        self.db_path = run_config.get("db_path", "")
        self.threads = run_config.get("threads", 1)
        self.extra = run_config.get("extra", "")

        from samovar.paths import resolve_executable

        raw_cmd = run_config.get("cmd") or self.default_cmd
        parts = command_argv(raw_cmd)
        # Prepared runs may retain the interpreter from an older installation.
        if len(parts) >= 3 and parts[1] == "-m" and parts[2].startswith("samovar."):
            raw_cmd = shlex.join([python_path(), *parts[1:]])
        first = command_argv(raw_cmd)[0] if raw_cmd else ""
        self.cmd = shlex.join(command_argv(resolve_executable(raw_cmd, tool_key=Path(first).name if first else None)))

    @property
    def default_cmd(self) -> str:
        raise NotImplementedError(
            "Method default_cmd must be implemented in child class"
        )

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        raise NotImplementedError(
            "Method get_expected_outputs must be implemented in child class"
        )

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        raise NotImplementedError(
            "Method get_snakemake_shell_cmd must be implemented in child class"
        )

    def parse_output(self, file_path: str) -> pd.DataFrame:
        raise NotImplementedError(
            "Method parse_output must be implemented in child class"
        )


@lru_cache(maxsize=32)
def _metaphlan_modern_cli(command):
    import subprocess
    result = subprocess.run(command_argv(command) + ["--help"], check=True, capture_output=True, text=True)
    return "--db_dir" in result.stdout


class ProfileAnnotator(BaseAnnotator):
    """Sample-level contract, sharing tool lookup/resources with annotators."""

    @property
    def default_cmd(self):
        return self.run_config["type"]

    def get_expected_outputs(self, sample, output_dir):
        out = os.path.join(output_dir, f"{sample}_{self.run_name}.profile.csv")
        return [out, out + ".profile.json"]

    def get_snakemake_shell_cmd(self, input_r1, input_r2, outputs):
        out = outputs[0]
        native = Path(command_argv(self.cmd)[-1]).name.lower().startswith("metaphlan")
        native = native or str(self.run_config["type"]).lower() in {"metaphlan", "metaphlan4", "mpa", "mp4"}
        sample = Path(out).name[:-len(f"_{self.run_name}.profile.csv")]
        py = shlex.quote(python_path())
        ingest = (f"{py} -m samovar.profiles ingest --output {shlex.quote(out)} "
                  f"--sample {shlex.quote(sample)} --identity {shlex.quote(self.run_name)} "
                  f"--database {shlex.quote(self.db_path)}")
        if str(self.run_config['type']).lower() == 'phanta' or Path(command_argv(self.cmd)[-1]).name == 'run_phanta.py':
            run = (f"{py} -m samovar.phanta -i {shlex.quote(str(input_r1))} {_r2_flag('-I', input_r2)}"
                   f"-d {shlex.quote(self.db_path)} -o {shlex.quote(out + '.phanta.tsv')} -t {int(self.threads)} "
                   f"--sample {shlex.quote(sample)} --exec-path {shlex.quote(self.cmd)} "
                   f"--native-flags={shlex.quote(self.extra or '')}")
            return run + f" && {ingest} --source {shlex.quote(out + '.phanta.tsv')} --format phanta"
        if str(self.run_config["type"]).lower() in {"motus", "motus3"} or Path(command_argv(self.cmd)[-1]).name == "motus":
            run = (f"{py} -m samovar.profiles motus --source {shlex.quote(out + '.raw')} "
                   f"--input-r1 {shlex.quote(str(input_r1))} {_r2_flag('--input-r2', input_r2)}"
                   f"--database {shlex.quote(self.db_path)} --sample {shlex.quote(sample)} "
                   f"--executable {shlex.quote(self.cmd)} --threads {int(self.threads)} "
                   f"--native-flags={shlex.quote(self.extra or '')}")
            return run + f" && {ingest} --source {shlex.quote(out + '.raw')} --format motus"
        if str(self.run_config["type"]).lower() == "bracken" or Path(command_argv(self.cmd)[-1]).name == "bracken":
            import argparse
            options = argparse.ArgumentParser(add_help=False)
            options.add_argument("-r", type=int)
            options.add_argument("-l", choices=["S"], default="S")
            options.add_argument("-t", type=int, default=10)
            opts = options.parse_args(shlex.split(self.extra or ""))
            lengths = [int(p.name[len("database"):-len("mers.kmer_distrib")])
                       for p in Path(self.db_path).glob("database*mers.kmer_distrib")
                       if p.name[len("database"):-len("mers.kmer_distrib")].isdigit()]
            length = opts.r if opts.r is not None else (lengths[0] if len(lengths) == 1 else 100)
            if length < 1 or opts.t < 0:
                raise ValueError("Bracken needs positive -r read length and nonnegative -t threshold")
            if not (Path(self.db_path) / f"database{length}mers.kmer_distrib").is_file():
                raise FileNotFoundError(f"Bracken distribution for -r {length} missing; build --type bracken or select -r")
            kraken = Kraken2Annotator({**self.run_config, "type": "kraken2", "cmd": "kraken2", "extra": ""}, self.config)
            reports = kraken.get_expected_outputs(sample, str(Path(out).parent))
            command = kraken.get_snakemake_shell_cmd(input_r1, input_r2, reports)
            native_profile = ingest.replace("samovar.profiles ingest", "samovar.profiles bracken")
            return command + (f" && {native_profile} --source {shlex.quote(out + '.bracken.tsv')} "
                              f"--kraken-report {shlex.quote(reports[0])} --read-length {length} --threshold {opts.t}")
        if str(self.run_config["type"]).lower() == "sylph" or Path(command_argv(self.cmd)[-1]).name.lower() == "sylph":
            db = Path(self.db_path)
            if db.is_dir():
                indexes = sorted([*db.glob("*.syldb"), *db.glob("*.syl2db")])
                if len(indexes) != 1:
                    raise ValueError("Sylph needs exactly one sketch database, or select its file explicitly")
                db = indexes[0]
            taxonomy = Path(str(db) + ".taxonomy.tsv")
            if not db.is_file() or not taxonomy.is_file():
                raise FileNotFoundError(f"Sylph needs a sketch and {taxonomy}; use samovar build --type sylph")
            reads = (f"-1 {shlex.quote(str(input_r1))} -2 {shlex.quote(str(input_r2))}" if input_r2
                     else shlex.quote(str(input_r1)))
            raw = out + ".raw"
            extra = shlex.split(self.extra or "")
            manifest = db.parent / "sylph.manifest.json"
            if db.name == "genomes.syldb" and manifest.is_file() and not any(p == "-c" or p.startswith("-c=") for p in extra):
                import json
                extra += ["-c", str(json.loads(manifest.read_text())["compression"])]
            run = f"{self.cmd} profile {shlex.quote(str(db))} {reads} -t {int(self.threads)} -o {shlex.quote(raw)} {shlex.join(extra)}"
            return skip_empty_reads_cmd(input_r1, [raw], run) + (
                f" && {ingest} --source {shlex.quote(raw)} --format sylph --genome-taxonomy {shlex.quote(str(taxonomy))}")
        if native:
            db = Path(self.db_path)
            modern = _metaphlan_modern_cli(self.cmd)
            db_flag, map_flag = ("--db_dir", "--mapout") if modern else ("--bowtie2db", "--bowtie2out")
            extra_parts = shlex.split(self.extra or "")
            extra_parts = [db_flag if p == "--bowtie2db" else p for p in extra_parts]
            if db_flag in extra_parts:
                at = extra_parts.index(db_flag)
                if at + 1 >= len(extra_parts) or Path(extra_parts[at + 1]) != db:
                    raise ValueError("MetaPhlAn database flags disagree with selected database")
                del extra_parts[at:at + 2]
            extra = shlex.join(extra_parts)
            if not any(p == "--bowtie2_exe" or p.startswith("--bowtie2_exe=") for p in extra_parts):
                from samovar.paths import resolve_executable
                extra += " --bowtie2_exe " + shlex.quote(resolve_executable("bowtie2", tool_key="bowtie2"))
            if "--index" not in shlex.split(extra):
                prefixes = sorted({p.name.split(".1.bt2")[0] for p in db.glob("*.1.bt2*") if ".rev.1.bt2" not in p.name})
                if len(prefixes) != 1:
                    raise ValueError("MetaPhlAn profiler needs a pinned --index or exactly one installed Bowtie2 index")
                extra += " --index " + shlex.quote(prefixes[0])
            raw = out + ".raw"
            reads = ",".join(p for p in (str(input_r1), str(input_r2 or "")) if p)
            run = (f"{self.cmd} {shlex.quote(reads)} --input_type fastq --nproc {int(self.threads)} "
                   f"{db_flag} {shlex.quote(self.db_path)} {map_flag} {shlex.quote(out + '.bowtie2.bz2')} "
                   f"-o {shlex.quote(raw)} {extra}")
            return skip_empty_reads_cmd(input_r1, [raw], run) + f" && {ingest} --source {shlex.quote(raw)} --format metaphlan"
        from samovar.paths import resolve_executable
        argv = command_argv(resolve_executable(self.cmd, allow_python_script=True))
        if Path(argv[0]).suffix == ".py":
            argv = [python_path(), *argv]
        run = (f"{shlex.join(argv)} -i {shlex.quote(str(input_r1))} {_r2_flag('-I', input_r2)}"
               f"-d {shlex.quote(self.db_path)} -o {shlex.quote(out)} -t {int(self.threads)} {self.extra or ''}")
        return run + f" && {ingest} --source {shlex.quote(out)} --format abundance"


class KMCPAnnotator(BaseAnnotator):
    """Read-level native KMCP search → profile -B through the shared contract."""

    @property
    def default_cmd(self):
        return "kmcp"

    def get_expected_outputs(self, sample, output_dir):
        # Reuse the universal annotation/ingestion rule for seq + taxID tools.
        return [os.path.join(output_dir, f"{sample}_{self.run_name}.custom_kmcp.out")]

    def get_snakemake_shell_cmd(self, input_r1, input_r2, outputs):
        return (f"{shlex.quote(python_path())} -m samovar.kmcp "
                f"-i {shlex.quote(str(input_r1))} {_r2_flag('-I', input_r2)}"
                f"-d {shlex.quote(str(self.db_path))} -o {shlex.quote(outputs[0])} "
                f"-t {int(self.threads)} --exec-path {shlex.quote(self.cmd)} {self.extra or ''}")

    def parse_output(self, file_path):
        from samovar.parse_annotators import read_custom_raw
        return read_custom_raw(file_path)


class Kraken2Annotator(BaseAnnotator):
    """Annotator class for Kraken2 taxonomic classifier."""

    @property
    def default_cmd(self) -> str:
        return "kraken2"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        report_file = os.path.join(
            output_dir, f"{sample}_{self.run_name}.kraken2.report"
        )
        out_file = os.path.join(output_dir, f"{sample}_{self.run_name}.kraken2.out")
        return [report_file, out_file]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        report_file, out_file = outputs
        extra = str(self.extra or "")
        mapping = "" if "--memory-mapping" in extra else "--memory-mapping "
        common = (
            f"{self.cmd} "
            f"--use-names "
            f"{mapping}"
            f"--db {shlex.quote(str(self.db_path))} "
            f"--threads {self.threads} "
            f"--report {shlex.quote(str(report_file))} "
            f"--output {shlex.quote(str(out_file))} "
            f"{extra}"
        )
        return single_or_paired_reads_cmd(
            input_r1,
            input_r2,
            outputs,
            single_cmd=f"{common} {shlex.quote(str(input_r1))}",
            paired_cmd=f"{common} --paired {shlex.quote(str(input_r1))} {shlex.quote(str(input_r2))}",
        )

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = _read_table_or_empty(file_path)
        if df is None:
            return _empty_taxid_frame()
        df.columns = ["classified", "seq", "taxa", "length", "k-mer"]

        def extract_taxid(taxa_string):
            if pd.isna(taxa_string):
                return "0"
            match = re.search(r"(?<=taxid )[0-9]*", str(taxa_string))
            return match.group(0) if match else "0"

        df["taxID"] = df["taxa"].apply(extract_taxid)
        if "length" in df.columns:
            df["length"] = df["length"].astype(str).str.replace(r"\|.*", "", regex=True)
            return df[["seq", "taxID", "length"]]
        return df[["seq", "taxID"]]


# One kaiju-multi process once a cohort is larger than a single sample.
KAIJU_MULTI_MIN_SAMPLES = 2


def use_kaiju_multi(n_samples: int) -> bool:
    return int(n_samples or 0) >= KAIJU_MULTI_MIN_SAMPLES


def _comma_paths(paths: Sequence[str]) -> str:
    return ",".join(shlex.quote(str(path)) for path in paths)


def _expand_index(path: str) -> str:
    """Turn ``dir/*.fmi`` into the one index file the shell used to expand."""
    text = str(path)
    if "*" not in text and "?" not in text:
        return text
    hits = sorted(glob.glob(text))
    if len(hits) == 1:
        return hits[0]
    return text


def kaiju_multi_script(
    *,
    binary: str,
    nodes: str,
    database: str,
    r1: Sequence[str],
    r2: Sequence[str],
    outputs: Sequence[str],
    threads: int,
    extra: str = "",
    empty: Optional[Sequence[str]] = None,
) -> str:
    """Shell for one ``kaiju-multi`` run. Empty FASTQs are touched and left out."""
    lines = ["set -euo pipefail"]
    for path in outputs:
        parent = str(path).rsplit("/", 1)[0] if "/" in str(path) else "."
        lines.append(f"mkdir -p {shlex.quote(parent)}")
    for path in empty or []:
        lines.append(f"touch {shlex.quote(str(path))}")
    if not r1:
        return "\n".join(lines)
    paired = [str(path) for path in r2 if str(path or "").strip()]
    if paired and len(paired) != len(r1):
        raise ValueError("kaiju-multi needs every sample paired, or none of them")
    cmd = [
        shlex.quote(binary),
        "-t",
        shlex.quote(nodes),
        "-f",
        shlex.quote(_expand_index(database)),
        "-i",
        _comma_paths(r1),
        "-o",
        _comma_paths(outputs),
        "-z",
        str(max(1, int(threads))),
    ]
    if paired:
        cmd.extend(["-j", _comma_paths(paired)])
    tail = str(extra or "").strip()
    if tail:
        cmd.append(tail)
    lines.append(" ".join(cmd))
    return "\n".join(lines)


def wrap_shell(script: str) -> str:
    """One argv so Snakemake's ``2> log`` covers the whole script."""
    return "bash -euo pipefail -c " + shlex.quote(script)


class KaijuAnnotator(BaseAnnotator):
    """Annotator class for Kaiju taxonomic classifier."""

    @property
    def default_cmd(self) -> str:
        return "kaiju"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        out_file = os.path.join(output_dir, f"{sample}_{self.run_name}.kaiju.out")
        return [out_file]

    def _database_and_nodes(self) -> tuple:
        if self.run_config.get("db_name"):
            db_file_path = os.path.join(self.db_path, self.run_config["db_name"])
        elif str(self.db_path).endswith(".fmi"):
            db_file_path = self.db_path
        else:
            db_file_path = os.path.join(self.db_path, "*.fmi")
        # Kaiju keeps nodes.dmp beside the index dir (…/kaijudb/nodes.dmp,
        # index in …/kaijudb/refseq/), not inside the refseq folder.
        explicit = str(self.run_config.get("db_nodes") or "").strip()
        folder = Path(str(db_file_path).replace("*", "")).parent
        db_nodes = explicit or str(folder / "nodes.dmp")
        try:
            nodes_ok = Path(str(db_nodes)).is_file()
        except OSError:
            nodes_ok = False
        if not nodes_ok:
            parent_nodes = folder.parent / "nodes.dmp"
            try:
                if parent_nodes.is_file():
                    db_nodes = str(parent_nodes)
                    nodes_ok = True
            except OSError:
                nodes_ok = False
        if not nodes_ok:
            try:
                from samovar.taxdump import nodes_dmp

                found = nodes_dmp()
                if found is not None:
                    db_nodes = str(found)
            except Exception:
                pass
        return _expand_index(str(db_file_path)), str(db_nodes)

    def batch_shell(self, samples: Sequence[str], r1_dir: str, r2_dir: str, output_dir: str) -> str:
        """One ``kaiju-multi`` for every sample of this run."""
        from samovar.paths import resolve_executable
        from samovar.seqio import fastq_has_reads, find_fastq_mate

        database, nodes = self._database_and_nodes()
        binary = command_argv(resolve_executable("kaiju-multi") or "kaiju-multi")[0]
        active_r1: List[str] = []
        active_r2: List[str] = []
        active_out: List[str] = []
        empty: List[str] = []
        paired = True
        for sample in samples:
            r1 = find_fastq_mate(r1_dir, sample, "R1")
            r2 = find_fastq_mate(r2_dir, sample, "R2")
            out = os.path.join(output_dir, f"{sample}_{self.run_name}.kaiju.out")
            if r1 is None or not fastq_has_reads(str(r1)):
                empty.append(out)
                continue
            if r2 is None:
                paired = False
            active_r1.append(str(r1))
            active_r2.append(str(r2) if r2 is not None else "")
            active_out.append(out)
        r2_arg: List[str] = []
        if paired and active_r2 and all(active_r2):
            r2_arg = active_r2
        elif any(active_r2):
            lines = ["set -euo pipefail"]
            for path in empty:
                lines.append(f"touch {shlex.quote(str(path))}")
            for r1, r2, out in zip(active_r1, active_r2, active_out):
                lines.append(
                    self.get_snakemake_shell_cmd(r1, r2, [out])
                )
            return wrap_shell("\n".join(lines))
        script = kaiju_multi_script(
            binary=binary,
            nodes=nodes,
            database=database,
            r1=active_r1,
            r2=r2_arg,
            outputs=active_out,
            threads=int(self.threads),
            extra=self.extra or "",
            empty=empty,
        )
        return wrap_shell(script)

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        db_file_path, db_nodes = self._database_and_nodes()

        common = (
            f"{self.cmd} "
            f"-t {shlex.quote(str(db_nodes))} "
            f"-f {shlex.quote(str(db_file_path))} "
            f"-i {shlex.quote(str(input_r1))} "
            f"-z {self.threads} "
            f"-o {shlex.quote(str(out_file))} "
            f"{self.extra}"
        )
        return single_or_paired_reads_cmd(
            input_r1,
            input_r2,
            outputs,
            single_cmd=common,
            paired_cmd=f"{common} -j {shlex.quote(str(input_r2))}",
        )

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = _read_table_or_empty(file_path)
        if df is None:
            return _empty_taxid_frame()
        df.columns = ["classified", "seq", "taxID"]
        return df[["seq", "taxID"]]


class MetaPhlanAnnotator(BaseAnnotator):
    """Annotator class for the MetaPhlAn taxonomic classifier."""

    @property
    def default_cmd(self) -> str:
        return "metaphlan"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        out_file = os.path.join(output_dir, f"{sample}_{self.run_name}.metaphlan4.out")
        raw_file = os.path.join(output_dir, f"{sample}_{self.run_name}.metaphlan4.raw")
        bowtie_file = os.path.join(output_dir, f"{sample}_{self.run_name}.bowtie2.bz2")
        return [out_file, raw_file, bowtie_file]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file, raw_file, bowtie_file = outputs
        cmd = (
            f"{self.cmd} "
            f"--input_type fastq "
            f"--nproc {self.threads} "
            + (
                f"-1 {shlex.quote(str(input_r1))} -2 {shlex.quote(str(input_r2))} "
                if _has_r2(input_r2)
                else f"{shlex.quote(str(input_r1))} "
            )
            + f"--bowtie2out {shlex.quote(str(bowtie_file))} "
            f"-o {shlex.quote(str(raw_file))} "
            f"{self.extra} "
            f"&& "
            f"bzcat {shlex.quote(str(bowtie_file))} | awk -F'\\t' '{{print $1 \"\\t\" $3}}' > {shlex.quote(str(out_file))}"
        )
        return skip_empty_reads_cmd(input_r1, outputs, cmd)

    def _get_db_mapping(self) -> Dict[str, str]:
        if not self.db_path:
            return {}
        from samovar.parse_annotators import resolve_metaphlan_db_file

        try:
            db_file = resolve_metaphlan_db_file(
                self.db_path, db_name=self.run_config.get("db_name")
            )
        except FileNotFoundError:
            return {}
        if not os.path.exists(db_file):
            return {}
        conn = sqlite3.connect(db_file)
        cursor = conn.cursor()
        cursor.execute("SELECT ref_id, tax_id FROM mpa_species_map")
        mapping = {str(row[0]): str(row[1]) for row in cursor.fetchall()}
        conn.close()
        return mapping

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = pd.read_table(file_path, header=None)
        df.columns = ["seq", "taxID"]
        db_mapping = self._get_db_mapping()

        if db_mapping:

            def extract_ref_id(tax_string: str) -> Optional[str]:
                if pd.isna(tax_string):
                    return None
                match = re.search(r"M\d+-c\d+", str(tax_string))
                return match.group(0) if match else None

            df["ref_id"] = df["taxID"].apply(extract_ref_id)
            df["taxID"] = df["ref_id"].map(db_mapping).fillna("0")
            df = df.drop(columns=["ref_id"])

        return df[["seq", "taxID"]]


class Kraken1Annotator(BaseAnnotator):
    """Annotator class for the original Kraken (Kraken 1) classifier."""

    @property
    def default_cmd(self) -> str:
        return "kraken"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        out_file = os.path.join(output_dir, f"{sample}_{self.run_name}.kraken.out")
        return [out_file]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        reads = (
            f"--paired {shlex.quote(str(input_r1))} {shlex.quote(str(input_r2))}"
            if _has_r2(input_r2)
            else shlex.quote(str(input_r1))
        )
        cmd = (
            f"{self.cmd} "
            f"--db {shlex.quote(str(self.db_path))} "
            f"--threads {self.threads} "
            f"{reads} "
            f"--output {shlex.quote(str(out_file))} "
            f"{self.extra}"
        )
        return skip_empty_reads_cmd(input_r1, outputs, cmd)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = _read_table_or_empty(file_path)
        if df is None:
            return _empty_taxid_frame()
        df.columns = ["classified", "seq", "taxID", "length", "k-mer"]
        df["length"] = df["length"].astype(str).str.replace(r"\|.*", "", regex=True)
        return df[["seq", "taxID", "length"]]


class KrakenUniqAnnotator(BaseAnnotator):
    """Annotator class for the KrakenUniq taxonomic classifier."""

    @property
    def default_cmd(self) -> str:
        return "krakenuniq"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        out_file = os.path.join(output_dir, f"{sample}_{self.run_name}.krakenuniq.out")
        return [out_file]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        reads = (
            f"--paired {shlex.quote(str(input_r1))} {shlex.quote(str(input_r2))}"
            if _has_r2(input_r2)
            else shlex.quote(str(input_r1))
        )
        cmd = (
            f"{self.cmd} "
            f"--db {shlex.quote(str(self.db_path))} "
            f"--threads {self.threads} "
            f"{reads} "
            f"--output {shlex.quote(str(out_file))} "
            f"{self.extra}"
        )
        return skip_empty_reads_cmd(input_r1, outputs, cmd)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = _read_table_or_empty(file_path)
        if df is None:
            return _empty_taxid_frame()
        df.columns = ["classified", "seq", "taxID", "length", "k-mer"]
        df["length"] = df["length"].astype(str).str.replace(r"\|.*", "", regex=True)
        return df[["seq", "taxID", "length"]]


def _cmd_is_custom_wrapper(cmd: str) -> bool:
    """True if ``cmd`` already points at custom.sh or a tool-specific wrapper."""
    if not cmd:
        return False
    text = str(cmd)
    if "custom.sh" in text:
        return True
    token = command_argv(text)[-1]
    return os.path.basename(token) in CUSTOM_WRAPPER_SCRIPTS


class CustomAnnotator(BaseAnnotator):
    """Generic Annotator class for any tool handled by custom.sh."""

    def __init__(self, run_config: Dict, config: Dict, tool_name: str):
        super().__init__(run_config, config)
        self.tool_name = tool_name

    @property
    def default_cmd(self) -> str:
        script = Path(__file__).resolve().parent.parent / "annotators" / "custom.sh"
        return f"bash {shlex.quote(str(script))}"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        out_file = os.path.join(
            output_dir, f"{sample}_{self.run_name}.custom_{self.tool_name}.out"
        )
        return [out_file]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        extra = self.extra or ""
        from samovar.tool_spec import threads_cli_token
        from samovar.main_config import lookup_tool_record
        from samovar.paths import load_config

        rec = lookup_tool_record(load_config(), self.tool_name)
        tflag = threads_cli_token(self.tool_name, self.threads, rec)
        from samovar.memory import translated_mem_flag

        raw_limit = self.config.get("mem_mb")
        mflag = translated_mem_flag(
            self.tool_name,
            self.threads,
            cores=int(self.config.get("cores") or 0) or None,
            limit_mb=int(raw_limit) if raw_limit else None,
            record=rec,
            existing=extra,
        )
        use_router = _cmd_is_custom_wrapper(self.cmd) or (
            str(self.tool_name).lower() in CUSTOM_SH_ROUTER_TOOLS
        )
        if use_router:
            cmd = self.cmd if _cmd_is_custom_wrapper(self.cmd) else self.default_cmd
            return (
                f"{cmd} "
                f"-i {shlex.quote(str(input_r1))} "
                f"{_r2_flag('-I', input_r2)}"
                f"-d {shlex.quote(str(self.db_path))} "
                f"-o {shlex.quote(str(out_file))} "
                f"-p {self.tool_name} "
                f"{tflag} "
                f"{mflag} "
                f"{extra}"
            )
        first = command_argv(self.cmd)[0]
        try:
            direct = bool(first) and Path(first).expanduser().is_file()
        except OSError:
            direct = False
        if direct:
            return (
                f"{self.cmd} "
                f"-i {shlex.quote(str(input_r1))} "
                f"{_r2_flag('-I', input_r2)}"
                f"-d {shlex.quote(str(self.db_path))} "
                f"-o {shlex.quote(str(out_file))} "
                f"{tflag} "
                f"{mflag} "
                f"{extra}"
            )
        cmd = self.default_cmd
        return (
            f"{cmd} "
            f"-i {shlex.quote(str(input_r1))} "
            f"{_r2_flag('-I', input_r2)}"
            f"-d {shlex.quote(str(self.db_path))} "
            f"-o {shlex.quote(str(out_file))} "
            f"-p {self.tool_name} "
            f"{tflag} "
            f"{mflag} "
            f"{extra}"
        )

    def parse_output(self, file_path: str) -> pd.DataFrame:
        from samovar.parse_annotators import read_custom_raw

        return read_custom_raw(file_path)


class ConstantTaxidAnnotator(CustomAnnotator):
    """Built-in classifier that assigns one taxID to every sequence.

    Default taxID is 9606 (Homo sapiens). Output uses the custom two-column
    seq/taxID table so it is consumed by the same Snakemake `custom_tool` rule.
    """

    def __init__(self, run_config: Dict, config: Dict, tool_name: str = "constant_taxid"):
        super().__init__(run_config, config, tool_name=tool_name)
        self.taxid = str(run_config.get("taxid", DUMMY_TAXID))

    @property
    def default_cmd(self) -> str:
        script = Path(__file__).resolve().parent / "baselines" / "constant_taxid.py"
        return f"{shlex.quote(python_path())} {shlex.quote(str(script))}"

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        cmd = self.cmd
        if not cmd or os.path.basename(command_argv(cmd)[0]).split(".")[0].lower() in CONSTANT_TAXID_TOOL_NAMES:
            cmd = self.default_cmd
        extra = self.extra or ""
        run = (
            f"{cmd} "
            f"-i {shlex.quote(str(input_r1))} "
            f"{_r2_flag('-I', input_r2)}"
            f"-o {shlex.quote(str(out_file))} "
            f"--taxid {shlex.quote(self.taxid)} "
            f"{extra}"
        )
        return skip_empty_reads_cmd(input_r1, [out_file], run)


class AssemblyAnnotator(BaseAnnotator):
    """FASTQ → assembly / MAG taxonomy / read taxIDs (nested contracts)."""

    @property
    def default_cmd(self) -> str:
        return f"{shlex.quote(python_path())} -m samovar.assembly_profiling"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        return [os.path.join(output_dir, f"{sample}_{self.run_name}.assembly.out")]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        sample_dir = os.path.join(
            os.path.dirname(out_file),
            os.path.basename(out_file).replace(".out", ""),
        )
        extra = self.extra or ""
        cmd = self.cmd or self.default_cmd
        # Prepare stores the bare name (`assembly_profiling`). That is not a
        # PATH binary; the implementation is this module.
        bare = os.path.basename(command_argv(cmd)[-1]) if str(cmd).strip() else ""
        if "samovar.assembly_profiling" not in str(cmd) and bare in {
            "assembly",
            "assembly_profiling",
        }:
            cmd = self.default_cmd
        run = (
            f"{cmd} annotator "
            f"-i {shlex.quote(str(input_r1))} "
            f"{_r2_flag('-I', input_r2)}"
            f"-d {shlex.quote(str(self.db_path or ''))} "
            f"-o {shlex.quote(str(out_file))} "
            f"-t {int(self.threads)} "
            f"--work-dir {shlex.quote(sample_dir)} "
            f"{extra}"
        )
        return skip_empty_reads_cmd(input_r1, [out_file], run)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        from samovar.parse_annotators import read_custom_raw

        return read_custom_raw(file_path)


GC_TOOL_NAMES = {"gc", "gc_content", "gc-content", "gcfrac", "gc_fraction"}
LENGTH_TOOL_NAMES = {"length", "read_length", "readlen", "read-length"}


class GcAnnotator(BaseAnnotator):
    """Per-read GC fraction Feature extractor."""

    @property
    def default_cmd(self) -> str:
        return f"{shlex.quote(python_path())} -m samovar.gc"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        return [os.path.join(output_dir, f"{sample}_{self.run_name}.gc.out")]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        extra = self.extra or ""
        cmd = self.cmd
        first = os.path.basename(command_argv(cmd)[0]) if cmd else ""
        if not cmd or first.split(".")[0].lower().replace("-", "_") in GC_TOOL_NAMES:
            cmd = self.default_cmd
        run = (
            f"{cmd} "
            f"-i {shlex.quote(str(input_r1))} "
            f"{_r2_flag('-I', input_r2)}"
            f"-o {shlex.quote(str(out_file))} "
            f"-t {int(self.threads)} "
            f"{extra}"
        )
        return skip_empty_reads_cmd(input_r1, [out_file], run)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        from samovar.parse_annotators import read_custom_raw

        return read_custom_raw(file_path)


class LengthAnnotator(BaseAnnotator):
    """Per-read length Feature extractor (R1 + R2 bases)."""

    @property
    def default_cmd(self) -> str:
        return f"{shlex.quote(python_path())} -m samovar.length"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        return [os.path.join(output_dir, f"{sample}_{self.run_name}.length.out")]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        extra = self.extra or ""
        cmd = self.cmd
        first = os.path.basename(command_argv(cmd)[0]) if cmd else ""
        if not cmd or first.split(".")[0].lower().replace("-", "_") in LENGTH_TOOL_NAMES:
            cmd = self.default_cmd
        run = (
            f"{cmd} "
            f"-i {shlex.quote(str(input_r1))} "
            f"{_r2_flag('-I', input_r2)}"
            f"-o {shlex.quote(str(out_file))} "
            f"-t {int(self.threads)} "
            f"{extra}"
        )
        return skip_empty_reads_cmd(input_r1, [out_file], run)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        from samovar.parse_annotators import read_custom_raw

        return read_custom_raw(file_path)


KMER2_TOOL_NAMES = {"kmer2", "kmer_counter", "dinuc", "dinucleotide"}


class Kmer2Annotator(BaseAnnotator):
    """Dinucleotide (k=2) feature extractor: FASTQ → seq + 16 ACGT 2-mer counts."""

    @property
    def default_cmd(self) -> str:
        return f"{shlex.quote(python_path())} -m samovar.kmer2"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        return [os.path.join(output_dir, f"{sample}_{self.run_name}.kmer2.out")]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        extra = self.extra or ""
        cmd = self.cmd
        first = os.path.basename(command_argv(cmd)[0]) if cmd else ""
        if not cmd or first.split(".")[0].lower() in KMER2_TOOL_NAMES:
            cmd = self.default_cmd
        run = (
            f"{cmd} "
            f"-i {shlex.quote(str(input_r1))} "
            f"{_r2_flag('-I', input_r2)}"
            f"-o {shlex.quote(str(out_file))} "
            f"-t {int(self.threads)} "
            f"{extra}"
        )
        return skip_empty_reads_cmd(input_r1, [out_file], run)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        from samovar.parse_annotators import read_custom_raw

        return read_custom_raw(file_path)


KMER_ENCODER_TOOL_NAMES = {
    "kmer_encoder",
    "kmer-encoder",
    "kmerenc",
    "encoder",
    "kmerencoder",
}


class KmerEncoderAnnotator(BaseAnnotator):
    """k-mer encoder Feature extractor (metauto encoding stage)."""

    @property
    def default_cmd(self) -> str:
        return f"{shlex.quote(python_path())} -m samovar.kmer_encoder"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        return [os.path.join(output_dir, f"{sample}_{self.run_name}.kmer_encoder.out")]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        extra = self.extra or ""
        cmd = self.cmd
        first = os.path.basename(command_argv(cmd)[0]) if cmd else ""
        if not cmd or first.split(".")[0].lower().replace("-", "_") in KMER_ENCODER_TOOL_NAMES:
            cmd = self.default_cmd
        db = self.db_path or ""
        run = (
            f"{cmd} "
            f"-i {shlex.quote(str(input_r1))} "
            f"{_r2_flag('-I', input_r2)}"
            f"-d {shlex.quote(str(db))} "
            f"-o {shlex.quote(str(out_file))} "
            f"-t {int(self.threads)} "
            f"{extra}"
        )
        return skip_empty_reads_cmd(input_r1, [out_file], run)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        from samovar.parse_annotators import read_custom_raw

        return read_custom_raw(file_path)


def get_annotator_instance(
    tool_type: str, run_config: Dict, config: Dict
) -> BaseAnnotator:
    """Factory function to instantiate the correct annotator class."""
    tool = tool_type.lower()
    if run_config.get("output_level") == "profile":
        return ProfileAnnotator(run_config, config)

    native_tools = {
        "kraken2": Kraken2Annotator,
        "kmcp": KMCPAnnotator,
        "kaiju": KaijuAnnotator,
        "metaphlan": MetaPhlanAnnotator,
        "metaphlan4": MetaPhlanAnnotator,
        "mpa": MetaPhlanAnnotator,
        "mp4": MetaPhlanAnnotator,
        "kraken": Kraken1Annotator,
        "kraken1": Kraken1Annotator,
        "krakenuniq": KrakenUniqAnnotator,
        "krakenu": KrakenUniqAnnotator,
        "assembly": AssemblyAnnotator,
        "assembly_profiling": AssemblyAnnotator,
        "gc": GcAnnotator,
        "gc_content": GcAnnotator,
        "gc-content": GcAnnotator,
        "gcfrac": GcAnnotator,
        "gc_fraction": GcAnnotator,
        "length": LengthAnnotator,
        "read_length": LengthAnnotator,
        "readlen": LengthAnnotator,
        "read-length": LengthAnnotator,
        "kmer2": Kmer2Annotator,
        "kmer_counter": Kmer2Annotator,
        "dinuc": Kmer2Annotator,
        "dinucleotide": Kmer2Annotator,
        "kmer_encoder": KmerEncoderAnnotator,
        "kmer-encoder": KmerEncoderAnnotator,
        "kmerenc": KmerEncoderAnnotator,
        "encoder": KmerEncoderAnnotator,
        "kmerencoder": KmerEncoderAnnotator,
    }

    parts = command_argv(run_config.get("cmd") or "")
    if parts and Path(parts[-1]).name == "kmcp":
        return KMCPAnnotator(run_config, config)
    if tool in CONSTANT_TAXID_TOOL_NAMES:
        return ConstantTaxidAnnotator(run_config, config, tool_name=tool_type)

    if tool in native_tools:
        return native_tools[tool](run_config, config)

    # Send any unknown tool to custom.sh
    return CustomAnnotator(run_config, config, tool_name=tool_type)
