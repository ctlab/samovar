import os
import re
import sqlite3
import sys
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

DUMMY_TAXID = "9606"
DUMMY_TOOL_NAMES = {
    "dummy",
    "dummy9606",
    "constant9606",
    "constant",
    "random",
}


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

        # Use custom command if provided in config, otherwise fall back to the tool's default
        self.cmd = run_config.get("cmd", self.default_cmd)

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
        cmd = (
            f"{self.cmd} "
            f"--use-names "
            f"--db {self.db_path} "
            f"--threads {self.threads} "
            f"--paired {input_r1} {input_r2} "
            f"--report {report_file} "
            f"--output {out_file} "
            f"{self.extra}"
        )
        return cmd

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = pd.read_table(file_path, header=None)
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
        db_nodes = self.run_config.get("db_nodes", default_nodes)

        cmd = (
            f"{self.cmd} "
            f"-t {db_nodes} "
            f"-f {db_file_path} "
            f"-i {input_r1} "
            f"-j {input_r2} "
            f"-z {self.threads} "
            f"-o {out_file} "
            f"{self.extra}"
        )
        return cmd

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = pd.read_table(file_path, header=None)
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
        return cmd

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
        return cmd

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = pd.read_table(file_path, header=None)
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
        return cmd

    def parse_output(self, file_path: str) -> pd.DataFrame:
        df = pd.read_table(file_path, header=None)
        df.columns = ["classified", "seq", "taxID", "length", "k-mer"]
        return df[["seq", "taxID"]]


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
        cmd = (
            f"{self.cmd} "
            f"-i {input_r1} "
            f"-I {input_r2} "
            f"-d {self.db_path} "
            f"-o {out_file} "
            f"-p {self.tool_name} "
            f"-t {self.threads} "
            f"{self.extra}"
        )
        return cmd

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
        return (
            f"{cmd} "
            f"-i {input_r1} "
            f"-I {input_r2} "
            f"-o {out_file} "
            f"--taxid {self.taxid} "
            f"{extra}"
        )


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
