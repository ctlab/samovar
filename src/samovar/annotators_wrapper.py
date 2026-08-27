import os
import re
import shlex
import sqlite3
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pandas as pd

DUMMY_TAXID = "9606"
DUMMY_TOOL_NAMES = {
    "dummy",
    "dummy9606",
    "constant9606",
    "constant",
    "random",
}

# Thin CLI wrappers that already append ``-p <tool>`` and forward custom.sh flags.
CUSTOM_WRAPPER_SCRIPTS = {
    "custom.sh",
    "centrifuge.sh",
    "metauto.sh",
    "assembly_hybrid.sh",
}

# These names always go through custom.sh (native CLIs are not the SamovaR flags).
CUSTOM_SH_ROUTER_TOOLS = {
    "centrifuge",
    "metauto",
    "assembly_hybrid",
} | DUMMY_TOOL_NAMES


def skip_empty_reads_cmd(input_r1: str, outputs: Sequence[str], run_cmd: str) -> str:
    """If R1 is empty, touch expected outputs instead of running the classifier.

    Kraken2/Kaiju on a 0-byte FASTQ can exit 0 without writing ``--output``.
    """
    r1 = shlex.quote(str(input_r1))
    outs = " ".join(shlex.quote(str(path)) for path in outputs)
    return f"if [ ! -s {r1} ]; then touch {outs}; else {run_cmd.strip()}; fi"


def single_or_paired_reads_cmd(
    input_r1: str,
    input_r2: str,
    outputs: Sequence[str],
    *,
    single_cmd: str,
    paired_cmd: str,
) -> str:
    """Run a classifier in paired mode only when R2 contains reads."""
    r1 = shlex.quote(str(input_r1))
    r2 = shlex.quote(str(input_r2))
    outs = " ".join(shlex.quote(str(path)) for path in outputs)
    return (
        f"if [ ! -s {r1} ]; then touch {outs}; "
        f"elif [ -s {r2} ]; then {paired_cmd.strip()}; "
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
        first = str(raw_cmd).split()[0] if raw_cmd else ""
        self.cmd = resolve_executable(raw_cmd, tool_key=Path(first).name if first else None)

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
        common = (
            f"{self.cmd} "
            f"--use-names "
            f"--db {self.db_path} "
            f"--threads {self.threads} "
            f"--report {report_file} "
            f"--output {out_file} "
            f"{self.extra}"
        )
        return single_or_paired_reads_cmd(
            input_r1,
            input_r2,
            outputs,
            single_cmd=f"{common} {input_r1}",
            paired_cmd=f"{common} --paired {input_r1} {input_r2}",
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
        return df[["seq", "taxID"]]


class KaijuAnnotator(BaseAnnotator):
    """Annotator class for Kaiju taxonomic classifier."""

    @property
    def default_cmd(self) -> str:
        return "kaiju"

    def get_expected_outputs(self, sample: str, output_dir: str) -> List[str]:
        out_file = os.path.join(output_dir, f"{sample}_{self.run_name}.kaiju.out")
        return [out_file]

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        if self.run_config.get("db_name"):
            db_file_path = os.path.join(self.db_path, self.run_config["db_name"])
        elif self.db_path.endswith(".fmi"):
            db_file_path = self.db_path
        else:
            db_file_path = os.path.join(self.db_path, "*.fmi")

        default_nodes = os.path.join(os.path.dirname(db_file_path), "nodes.dmp")
        db_nodes = self.run_config.get("db_nodes") or default_nodes
        try:
            nodes_ok = Path(str(db_nodes)).is_file()
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

        common = (
            f"{self.cmd} "
            f"-t {db_nodes} "
            f"-f {db_file_path} "
            f"-i {input_r1} "
            f"-z {self.threads} "
            f"-o {out_file} "
            f"{self.extra}"
        )
        return single_or_paired_reads_cmd(
            input_r1,
            input_r2,
            outputs,
            single_cmd=common,
            paired_cmd=f"{common} -j {input_r2}",
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
            f"-1 {input_r1} "
            f"-2 {input_r2} "
            f"--bowtie2out {bowtie_file} "
            f"-o {raw_file} "
            f"{self.extra} "
            f"&& "
            f"bzcat {bowtie_file} | awk -F'\\t' '{{print $1 \"\\t\" $3}}' > {out_file}"
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
        cmd = (
            f"{self.cmd} "
            f"--db {self.db_path} "
            f"--threads {self.threads} "
            f"--paired {input_r1} {input_r2} "
            f"--output {out_file} "
            f"{self.extra}"
        )
        return skip_empty_reads_cmd(input_r1, outputs, cmd)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = _read_table_or_empty(file_path)
        if df is None:
            return _empty_taxid_frame()
        df.columns = ["classified", "seq", "taxID", "length", "k-mer"]
        return df[["seq", "taxID"]]


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
        cmd = (
            f"{self.cmd} "
            f"--db {self.db_path} "
            f"--threads {self.threads} "
            f"--paired {input_r1} {input_r2} "
            f"--output {out_file} "
            f"{self.extra}"
        )
        return skip_empty_reads_cmd(input_r1, outputs, cmd)

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = _read_table_or_empty(file_path)
        if df is None:
            return _empty_taxid_frame()
        df.columns = ["classified", "seq", "taxID", "length", "k-mer"]
        return df[["seq", "taxID"]]


def _cmd_is_custom_wrapper(cmd: str) -> bool:
    """True if ``cmd`` already points at custom.sh or a tool-specific wrapper."""
    if not cmd:
        return False
    text = str(cmd)
    if "custom.sh" in text:
        return True
    token = text.split()[-1]
    return os.path.basename(token) in CUSTOM_WRAPPER_SCRIPTS


class CustomAnnotator(BaseAnnotator):
    """Generic Annotator class for any tool handled by custom.sh."""

    def __init__(self, run_config: Dict, config: Dict, tool_name: str):
        super().__init__(run_config, config)
        self.tool_name = tool_name

    @property
    def default_cmd(self) -> str:
        script = Path(__file__).resolve().parent.parent / "annotators" / "custom.sh"
        return f"bash {script}"

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
        use_router = _cmd_is_custom_wrapper(self.cmd) or (
            str(self.tool_name).lower() in CUSTOM_SH_ROUTER_TOOLS
        )
        if use_router:
            cmd = self.cmd if _cmd_is_custom_wrapper(self.cmd) else self.default_cmd
            return (
                f"{cmd} "
                f"-i {input_r1} "
                f"-I {input_r2} "
                f"-d {self.db_path} "
                f"-o {out_file} "
                f"-p {self.tool_name} "
                f"-t {self.threads} "
                f"{extra}"
            )
        first = str(self.cmd or "").split()[0]
        try:
            direct = bool(first) and Path(first).expanduser().is_file()
        except OSError:
            direct = False
        if direct:
            return (
                f"{self.cmd} "
                f"-i {input_r1} "
                f"-I {input_r2} "
                f"-d {self.db_path} "
                f"-o {out_file} "
                f"-t {self.threads} "
                f"{extra}"
            )
        cmd = self.default_cmd
        return (
            f"{cmd} "
            f"-i {input_r1} "
            f"-I {input_r2} "
            f"-d {self.db_path} "
            f"-o {out_file} "
            f"-p {self.tool_name} "
            f"-t {self.threads} "
            f"{extra}"
        )

    def parse_output(self, file_path: str) -> pd.DataFrame:
        # Safe reading for empty files
        try:
            if os.stat(file_path).st_size == 0:
                return pd.DataFrame(columns=["seq", "taxID"])
            df = pd.read_table(file_path, header=None).iloc[:, [0, 1]]
            df.columns = ["seq", "taxID"]
            return df
        except (pd.errors.EmptyDataError, FileNotFoundError):
            return pd.DataFrame(columns=["seq", "taxID"])


class ConstantTaxidAnnotator(CustomAnnotator):
    """Dummy / custom classifier that assigns one taxID to every sequence.

    Default taxID is 9606 (Homo sapiens). Output uses the custom two-column
    seq/taxID table so it is consumed by the same Snakemake `custom_tool` rule.
    """

    def __init__(self, run_config: Dict, config: Dict, tool_name: str = "constant9606"):
        super().__init__(run_config, config, tool_name=tool_name)
        self.taxid = str(run_config.get("taxid", DUMMY_TAXID))

    @property
    def default_cmd(self) -> str:
        script = Path(__file__).resolve().parent.parent / "annotators" / "constant9606.py"
        return f"{sys.executable} {script}"

    def get_snakemake_shell_cmd(
        self, input_r1: str, input_r2: str, outputs: List[str]
    ) -> str:
        out_file = outputs[0]
        cmd = self.cmd
        if not cmd or os.path.basename(str(cmd).split()[0]).split(".")[0].lower() in DUMMY_TOOL_NAMES:
            cmd = self.default_cmd
        extra = self.extra or ""
        run = (
            f"{cmd} "
            f"-i {shlex.quote(str(input_r1))} "
            f"-I {shlex.quote(str(input_r2))} "
            f"-o {shlex.quote(str(out_file))} "
            f"--taxid {shlex.quote(self.taxid)} "
            f"{extra}"
        )
        return skip_empty_reads_cmd(input_r1, [out_file], run)


def get_annotator_instance(
    tool_type: str, run_config: Dict, config: Dict
) -> BaseAnnotator:
    """Factory function to instantiate the correct annotator class."""
    tool = tool_type.lower()

    native_tools = {
        "kraken2": Kraken2Annotator,
        "kaiju": KaijuAnnotator,
        "metaphlan": MetaPhlanAnnotator,
        "metaphlan4": MetaPhlanAnnotator,
        "mpa": MetaPhlanAnnotator,
        "mp4": MetaPhlanAnnotator,
        "kraken": Kraken1Annotator,
        "kraken1": Kraken1Annotator,
        "krakenuniq": KrakenUniqAnnotator,
        "krakenu": KrakenUniqAnnotator,
    }

    if tool in DUMMY_TOOL_NAMES:
        return ConstantTaxidAnnotator(run_config, config, tool_name=tool_type)

    if tool in native_tools:
        return native_tools[tool](run_config, config)

    # Send any unknown tool to custom.sh
    return CustomAnnotator(run_config, config, tool_name=tool_type)
