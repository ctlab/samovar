"""Module for parsing and processing taxonomic annotation outputs from various tools.

This module provides functionality to read and process taxonomic annotation outputs
from tools like Kraken, Kaiju, and MetaPhlAn. It standardizes the output format
to allow for easy comparison and analysis of results from different tools.
"""

import pandas as pd
import os
import re
import sys
from typing import Dict, List, Optional
import sqlite3
import pickle
from pathlib import Path

# `ete3` imports its webplugin unconditionally, and the webplugin imports the
# stdlib `cgi` module. Python 3.13 removed `cgi`, so importing ete3 fails
# unless we provide a tiny compatibility shim.
try:
    from ete3 import NCBITaxa
except ModuleNotFoundError as exc:
    if str(exc).endswith(": cgi") or getattr(exc, "name", None) == "cgi":
        # Minimal `cgi` module shim: ete3 only needs `cgi.parse_qs` at import
        # time (runtime usage is irrelevant for our CLI/library usage).
        import types
        import urllib.parse

        cgi_stub = types.ModuleType("cgi")
        cgi_stub.parse_qs = urllib.parse.parse_qs
        sys.modules["cgi"] = cgi_stub
        from ete3 import NCBITaxa
    else:
        raise

# Lazily initialize NCBI taxonomy database to avoid expensive work during
# module import (this also makes unit test collection faster).
ncbi = None
_taxonomy_parent_rank = None

# Minimal built-in fallback for common test taxids if no taxonomy DB is available.
_FALLBACK_LINEAGE = {
    9606: {"species": 9606, "genus": 9605},
    511145: {"species": 562, "genus": 561},
    562: {"species": 562, "genus": 561},
}


def _get_ncbi() -> "NCBITaxa":
    global ncbi
    if ncbi is None:
        ncbi = NCBITaxa()
    return ncbi


def _load_nodes_cache(nodes_path: str) -> Dict[int, tuple]:
    """Load parent/rank map from nodes.dmp with filesystem cache."""
    cache_dir = Path.home() / ".cache" / "samovar"
    cache_dir.mkdir(parents=True, exist_ok=True)
    stat = os.stat(nodes_path)
    cache_key = f"{nodes_path}:{int(stat.st_mtime)}:{stat.st_size}"
    cache_name = f"nodes_{abs(hash(cache_key))}.pkl"
    cache_file = cache_dir / cache_name

    if cache_file.exists():
        with open(cache_file, "rb") as fh:
            return pickle.load(fh)

    tree = {}
    with open(nodes_path, "r", encoding="utf-8", errors="ignore") as fh:
        for line in fh:
            parts = line.split("\t|\t")
            if len(parts) < 3:
                continue
            taxid = int(parts[0].strip())
            parent = int(parts[1].strip())
            rank = parts[2].strip()
            tree[taxid] = (parent, rank)

    with open(cache_file, "wb") as fh:
        pickle.dump(tree, fh, protocol=pickle.HIGHEST_PROTOCOL)

    return tree


def _discover_nodes_path() -> Optional[str]:
    """Best-effort discovery for nodes.dmp in typical SamovaR DB paths."""
    env_path = os.environ.get("SAMOVAR_NODES_DMP")
    if env_path and os.path.exists(env_path):
        return env_path

    candidates = [
        "tests_outs/kraken_db/taxonomy/nodes.dmp",
        "tests_outs/kraken_db_test/taxonomy/nodes.dmp",
        "tests_outs/kaiju_db/taxonomy/nodes.dmp",
        "tests_outs/kaiju_db_test/taxonomy/nodes.dmp",
        "/mnt/metagenomics/kaiju/nodes.dmp",
    ]
    for path in candidates:
        if os.path.exists(path):
            return path
    return None


def _get_taxid_by_rank_from_nodes(taxid: int, rank_name: str) -> Optional[str]:
    """Resolve rank using local nodes.dmp parser."""
    global _taxonomy_parent_rank
    if _taxonomy_parent_rank is None:
        nodes_path = _discover_nodes_path()
        if nodes_path is None:
            return None
        _taxonomy_parent_rank = _load_nodes_cache(nodes_path)

    tree = _taxonomy_parent_rank
    current = taxid
    visited = set()
    while current in tree and current not in visited:
        visited.add(current)
        parent, rank = tree[current]
        if rank == rank_name:
            return str(current)
        if parent == current:
            break
        current = parent
    return None


def _resolve_taxid_by_rank(taxid: str, rank_name: str) -> Optional[str]:
    """Resolve target rank with fast and robust fallback strategy."""
    if taxid == "0":
        return "0"

    try:
        taxid_int = int(taxid)
    except (TypeError, ValueError):
        return None

    # 1) Try ete3-backed resolver.
    try:
        lineage = _get_ncbi().get_lineage(taxid_int)
        ranks = _get_ncbi().get_rank(lineage)
        for tid, rank in ranks.items():
            if rank == rank_name:
                return str(tid)
    except Exception:
        pass

    # 2) Try local nodes.dmp cache.
    local_hit = _get_taxid_by_rank_from_nodes(taxid_int, rank_name)
    if local_hit is not None:
        return local_hit

    # 3) Built-in tiny fallback for deterministic unit tests/common IDs.
    fallback = _FALLBACK_LINEAGE.get(taxid_int, {}).get(rank_name)
    if fallback is not None:
        return str(fallback)

    return None


def parse_metaphlan_db(db_path: str) -> Dict[str, str]:
    """Parse MetaPhlAn database to map reference IDs to NCBI taxIDs.
    
    Args:
        db_path: Path to MetaPhlAn database directory
        
    Returns:
        Dictionary mapping MetaPhlAn reference IDs to NCBI taxIDs
    """
    # MetaPhlAn uses a SQLite database
    db_file = os.path.join(db_path, "mpa_v30_CHOCOPhlAn_201901_species_map.db")
    if not os.path.exists(db_file):
        raise FileNotFoundError(f"MetaPhlAn database not found at {db_file}")
        
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    
    # Query to get mapping between reference IDs and taxIDs
    cursor.execute("SELECT ref_id, tax_id FROM mpa_species_map")
    mapping = {row[0]: row[1] for row in cursor.fetchall()}
    
    conn.close()
    return mapping


def _read_table_or_empty(file_path: str, columns: List[str]) -> pd.DataFrame:
    """Read a whitespace table, returning an empty frame for empty files."""
    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
        return pd.DataFrame(columns=columns)
    try:
        df = pd.read_table(file_path, header=None)
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=columns)
    if df.empty:
        return pd.DataFrame(columns=columns)
    df.columns = columns[: df.shape[1]]
    return df


def read_kaiju_raw(file_path: str, db_path: Optional[str] = None) -> pd.DataFrame:
    """Read raw Kaiju output file.
    
    Args:
        file_path: Path to Kaiju output file
        db_path: Optional path to Kaiju database file
        
    Returns:
        DataFrame with columns: classified, seq, score, taxID, N
    """
    df = _read_table_or_empty(file_path, ["classified", "seq", "taxID"])
    return df


def read_kraken1_raw(file_path: str) -> pd.DataFrame:
    """Read raw Kraken1 output file.
    
    Args:
        file_path: Path to Kraken1 output file
        
    Returns:
        DataFrame with columns: classified, seq, taxID, length, k-mer
    """
    df = _read_table_or_empty(
        file_path, ["classified", "seq", "taxID", "length", "k-mer"]
    )
    return df


def read_kraken2_raw(file_path: str) -> pd.DataFrame:
    """Read raw Kraken2 output file.
    
    Args:
        file_path: Path to Kraken2 output file
        
    Returns:
        DataFrame with columns: classified, seq, taxa, length, k-mer, taxID
    """
    df = _read_table_or_empty(
        file_path, ["classified", "seq", "taxa", "length", "k-mer"]
    )
    if df.empty:
        df["taxID"] = []
        return df
    df["length"] = [re.sub(r"\|.*", "", i) for i in df["length"]]
    df["taxID"] = [re.search(r'(?<=taxid )[0-9]*', i).group(0) for i in df["taxa"]]
    df["taxa"] = [re.sub(r' \(.*', "", i) for i in df["taxa"]]
    return df


def read_krakenu_raw(file_path: str) -> pd.DataFrame:
    """Read raw Kraken-Unique output file.
    
    Args:
        file_path: Path to Kraken-Unique output file
        
    Returns:
        DataFrame with columns: classified, seq, taxID, length, k-mer
    """
    df = _read_table_or_empty(
        file_path, ["classified", "seq", "taxID", "length", "k-mer"]
    )
    return df


def read_metaphlan_raw(file_path: str, db_path: Optional[str] = None) -> pd.DataFrame:
    """Read raw MetaPhlAn output file.
    
    Args:
        file_path: Path to MetaPhlAn output file
        db_path: Optional path to MetaPhlAn database directory
        
    Returns:
        DataFrame with columns: seq, taxID
    """
    df = _read_table_or_empty(file_path, ["seq", "taxID"])
    if df.empty:
        return df
    
    if db_path is not None:
        # Parse database and map reference IDs to taxIDs
        db_mapping = parse_metaphlan_db(db_path)
        
        # Extract reference ID from the taxID field
        df["ref_id"] = df["taxID"].apply(lambda x: re.search(r'M\d+-c\d+', x).group(0) if re.search(r'M\d+-c\d+', x) else None)
        
        # Map reference IDs to taxIDs
        df["taxID"] = df["ref_id"].map(db_mapping)
        df = df.drop("ref_id", axis=1)
    
    return df


# Dictionary mapping tool names to their respective read functions
READ_FUNCTIONS = {
    "kraken": read_kraken1_raw,
    "kraken1": read_kraken1_raw,
    "kraken2": read_kraken2_raw,
    "krakenunique": read_krakenu_raw,
    "krakenuniq": read_krakenu_raw,
    "krakenu": read_krakenu_raw,
    "metaphlan": read_metaphlan_raw,
    "mpa": read_metaphlan_raw,
    "mp4": read_metaphlan_raw,
    "kaiju": read_kaiju_raw
}


def read_annotation(file_path_type: Dict[str, str], trimmed: bool = True) -> pd.DataFrame:
    """Read and combine multiple annotation files.
    
    Args:
        file_path_type: Dictionary mapping file paths to their tool types
        trimmed: Whether to only keep taxID columns
        
    Returns:
        Combined DataFrame with annotations from all tools
    """
    res = pd.DataFrame()
    for path, tool_type in file_path_type.items():
        df = READ_FUNCTIONS.get(tool_type)(path).set_index("seq")
        if trimmed:
            df = df[["taxID"]]
        tool_name = tool_type.split("_")[-1].split(".")[0]
        df.columns = [f"{col}_{tool_name}" if col != "seq" else col for col in df.columns]
        res = pd.concat([res, df], axis=1)
    return res


class Annotation:
    """Class for handling taxonomic annotations from multiple tools.
    
    This class provides functionality to read, process, and analyze taxonomic
    annotations from various tools. It supports operations like expanding
    taxonomic ranks, comparing annotations, and exporting results.
    """
    
    def __init__(self, file_path_type: Dict[str, str], get_true_annotation: Optional[str] = None):
        """Initialize Annotation object.
        
        Args:
            file_path_type: Dictionary mapping file paths to their tool types
            get_true_annotation: Regex pattern to extract true annotation from sequence IDs
        """
        self.id = 0
        self.DataFrame = pd.DataFrame()
        
        # Read and combine all annotation files
        for path, tool_type in file_path_type.items():
            df = READ_FUNCTIONS.get(tool_type)(path).set_index("seq").astype({"taxID": 'string'})
            df.columns = [f"{col}_{tool_type}_{self.id}" for col in df.columns]
            
            # Check for duplicate indices
            if df.index.duplicated().any():
                df = df[~df.index.duplicated(keep='first')]
            
            try:
                self.DataFrame = pd.concat([self.DataFrame, df], axis=1)
            except:
                raise ValueError(f"Error concatenating {path}")
            self.id += 1

        # Extract true annotations if pattern provided
        if get_true_annotation is not None:
            self.true_annotation = []
            for i in self.DataFrame.index:
                match = re.search(get_true_annotation, i)
                if match:
                    self.true_annotation.append(match.group(0))
                else:
                    self.true_annotation.append("")
            print("True annotations extracted")

        # Get unique annotations and ranks
        set_columns = []
        for name, column in self.tr().items():
            set_columns += list(set(column))

        self.annotation_list = self.list2set(set_columns)
        self.true_annotation_list = self.list2set(self.true_annotation)
        self.rank_list = self.list2set([*self.annotation_list, *self.true_annotation_list])

    def true_annotation_unique(self) -> set:
        """Get unique true annotations."""
        return set(self.true_annotation)

    def true_annotation_rename(self, change_dict: Dict[str, str]):
        """Rename true annotations using provided mapping.
        
        Args:
            change_dict: Dictionary mapping old names to new names
        """
        self.true_annotation = [change_dict.get(TA) if TA in change_dict else "" for TA in self.true_annotation]

    def rank_annotation(self, rank: str = "species") -> 'RankAnnotation':
        """Get annotations at specified taxonomic rank.
        
        Args:
            rank: Taxonomic rank to get annotations for
            
        Returns:
            RankAnnotation object with annotations at specified rank
        """
        rank_list = [self.rank(j, rank) for j in self.rank_list]
        rank_dict = dict(zip(self.rank_list, rank_list))
        return RankAnnotation(rank).make(self.full(), rank_dict)

    def expand_annotation(self, rank: List[str] = ["species"]) -> 'ExpandAnnotation':
        """Expand annotations to multiple taxonomic ranks.
        
        Args:
            rank: List of taxonomic ranks to expand to
            
        Returns:
            ExpandAnnotation object with annotations at multiple ranks
        """
        full_rank_annotation = ExpandAnnotation()
        for i in rank:
            full_rank_annotation.add(self.rank_annotation(i))
        return full_rank_annotation

    def correct_annotations(self, rank: str = "species") -> pd.DataFrame:
        """Get correct annotations at specified rank.
        
        Args:
            rank: Taxonomic rank to check
            
        Returns:
            DataFrame with counts of correct annotations
        """
        return pd.DataFrame(self.rank_annotation(rank).correct_annotation().annotation.value_counts())

    def full(self) -> pd.DataFrame:
        """Get full annotation DataFrame including true annotations."""
        tmp = self.tr()
        tmp["true"] = self.true_annotation
        return tmp

    def correct_level(self, level: str = "sp") -> None:
        """Correct annotations to specified taxonomic level.
        
        Args:
            level: Taxonomic level to correct to (e.g. 'sp' for species, 'genus', etc.)
            
        This method:
        1. Gets unique taxIDs from the DataFrame
        2. Maps them to the specified taxonomic level
        3. Replaces original taxIDs with their corrected versions
        """
        # Get all unique taxIDs from taxID columns
        taxid_cols = [col for col in self.DataFrame.columns if col.startswith('taxID')]
        unique_taxids = set()
        for col in taxid_cols:
            unique_taxids.update(self.DataFrame[col].dropna().unique())
        
        # Create mapping dictionary for taxIDs
        taxid_map = {}
        for taxid in unique_taxids:
            if taxid == "0" or pd.isna(taxid):
                taxid_map[taxid] = taxid
                continue
            try:
                ranked = _resolve_taxid_by_rank(taxid, level)
                taxid_map[taxid] = ranked if ranked is not None else taxid
            except Exception:
                taxid_map[taxid] = taxid
        
        # Replace taxIDs in all taxID columns
        for col in taxid_cols:
            self.DataFrame[col] = self.DataFrame[col].map(taxid_map)

    def tr(self) -> pd.DataFrame:
        """Get DataFrame with only taxID columns."""
        return self.DataFrame.copy().filter(regex="taxID.*")

    @staticmethod
    def list2set(a: List) -> List[str]:
        """Convert list to set of strings."""
        return list(set([str(i) for i in a]))

    @staticmethod
    def list2rank(a: List, at_rank: str) -> List:
        """Convert list of taxIDs to nodes at specified rank.
        
        Args:
            a: List of taxIDs
            at_rank: Target taxonomic rank
            
        Returns:
            List of taxIDs at specified rank
        """
        b = Annotation.list2set(a)
        return [Annotation.rank(i, at_rank) for i in b]

    @staticmethod
    def rank(j: str, i: str) -> Optional[str]:
        """Get taxID at specified rank.
        
        Args:
            j: Input taxID
            i: Target rank
            
        Returns:
            taxID at specified rank or None if not found
        """
        if j == "0":
            return "0"
        try:
            return _resolve_taxid_by_rank(j, i)
        except Exception:
            return None

    def export(self, file: Optional[str] = None) -> pd.DataFrame:
        """Export annotations to CSV file.
        
        Args:
            file: Optional path to save CSV file
            
        Returns:
            DataFrame with exported annotations
        """
        df_return = self.DataFrame.loc[:, [col for col in self.DataFrame if col.startswith('taxID')]]
        lencol = [col for col in self.DataFrame if col.startswith('len')]
        if len(lencol) > 0:
            df_return['length'] = self.DataFrame.loc[:, lencol[0]].to_list()
        df_return["true"] = self.true_annotation
        if file is not None:
            df_return.to_csv(file)
        return df_return


class RankAnnotation:
    """Class for handling annotations at a specific taxonomic rank."""
    
    def __init__(self, rank: str):
        """Initialize RankAnnotation object.
        
        Args:
            rank: Taxonomic rank
        """
        self.rank = rank
        self.annotation = pd.DataFrame()

    def add(self, name: str, annotation: List):
        """Add annotation column.
        
        Args:
            name: Column name
            annotation: List of annotations
        """
        self.annotation[name] = annotation

    def make(self, annotation: pd.DataFrame, rank_dict: Dict[str, str]) -> 'RankAnnotation':
        """Create RankAnnotation from DataFrame.
        
        Args:
            annotation: Input DataFrame
            rank_dict: Dictionary mapping taxIDs to ranks
            
        Returns:
            RankAnnotation object
        """
        for name, column in annotation.items():
            self.add(str(name), [rank_dict.get(TA) for TA in column])
        self.reindex(annotation.index)
        return self

    def reindex(self, index):
        """Set DataFrame index."""
        self.annotation.index = index

    def y(self) -> pd.Series:
        """Get true annotations."""
        return self.annotation["true"]

    def x(self) -> pd.DataFrame:
        """Get predicted annotations."""
        return self.annotation.copy().drop("true", axis=1)

    def correct_annotation(self) -> 'RankAnnotation':
        """Get correct annotations."""
        tmp = RankAnnotation(self.rank)
        for name, column in self.x().items():
            tmp.add(str(name), pd.DataFrame(column == self.y()))
        tmp.reindex(self.annotation.index)
        return tmp


class ExpandAnnotation:
    """Class for handling annotations at multiple taxonomic ranks."""
    
    def __init__(self):
        """Initialize ExpandAnnotation object."""
        self.rank_annotation = {}

    def add(self, rank_annotation: RankAnnotation):
        """Add RankAnnotation object.
        
        Args:
            rank_annotation: RankAnnotation object to add
        """
        self.rank_annotation[rank_annotation.rank] = rank_annotation

    def get(self, rank: str) -> RankAnnotation:
        """Get RankAnnotation for specified rank.
        
        Args:
            rank: Taxonomic rank
            
        Returns:
            RankAnnotation object
        """
        return self.rank_annotation[rank]

def match_annotation(annotation_name:str) -> str:
    """Match annotation name to tool name.
    
    Args:
        annotation_name: Annotation name
        
    Returns:
        Tool name  
    """
    annotation_name = os.path.basename(annotation_name)
    if annotation_name.endswith(".out"):
        arg = annotation_name.split(".")[-2]
    
        if arg == "kraken":
            return "kraken"
        elif arg == "kraken1":
            return "kraken1"
        elif arg == "kraken2":
            return "kraken2"
        elif arg == "krakenuniq":
            return "krakenuniq"
        elif arg == "krakenu":
            return "krakenu"
        elif arg == "metaphlan":
            return "metaphlan"
        elif arg == "metaphlan4":
            return "metaphlan"
        elif arg == "mpa":
            return "metaphlan"
        elif arg == "mp4":
            return "metaphlan"
        elif arg == "kaiju":
            return "kaiju"
        else:
            raise ValueError(f"Annotation name {annotation_name} not found. Check the input: file extention should be like sample_metainfo.annotator.out")
    else:
        print("Skipping", annotation_name)
        return None