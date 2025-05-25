"""Module for parsing and processing taxonomic annotation outputs from various tools.

This module provides functionality to read and process taxonomic annotation outputs
from tools like Kraken, Kaiju, and MetaPhlAn. It standardizes the output format
to allow for easy comparison and analysis of results from different tools.
"""

import pandas as pd
import os
import re
from ete3 import NCBITaxa
from typing import Dict, List, Optional, Union
import sqlite3
import mmap

# Initialize NCBI taxonomy database
ncbi = NCBITaxa()


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


def read_kaiju_raw(file_path: str, db_path: Optional[str] = None) -> pd.DataFrame:
    """Read raw Kaiju output file.
    
    Args:
        file_path: Path to Kaiju output file
        db_path: Optional path to Kaiju database file
        
    Returns:
        DataFrame with columns: classified, seq, score, taxID, N
    """
    df = pd.read_table(file_path, header=None)
    df.columns = ["classified", "seq", "taxID"]
    return df


def read_kraken1_raw(file_path: str) -> pd.DataFrame:
    """Read raw Kraken1 output file.
    
    Args:
        file_path: Path to Kraken1 output file
        
    Returns:
        DataFrame with columns: classified, seq, taxID, length, k-mer
    """
    df = pd.read_table(file_path, header=None)
    df.columns = ["classified", "seq", "taxID", "length", "k-mer"]
    return df


def read_kraken2_raw(file_path: str) -> pd.DataFrame:
    """Read raw Kraken2 output file.
    
    Args:
        file_path: Path to Kraken2 output file
        
    Returns:
        DataFrame with columns: classified, seq, taxa, length, k-mer, taxID
    """
    df = pd.read_table(file_path, header=None)
    df.columns = ["classified", "seq", "taxa", "length", "k-mer"]
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
    df = pd.read_table(file_path, header=None)
    df.columns = ["classified", "seq", "taxID", "length", "k-mer"]
    return df


def read_metaphlan_raw(file_path: str, db_path: Optional[str] = None) -> pd.DataFrame:
    """Read raw MetaPhlAn output file.
    
    Args:
        file_path: Path to MetaPhlAn output file
        db_path: Optional path to MetaPhlAn database directory
        
    Returns:
        DataFrame with columns: seq, taxID
    """
    df = pd.read_table(file_path, header=None)
    df.columns = ["seq", "taxID"]
    
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
            try:
                self.DataFrame = pd.concat([self.DataFrame, df], axis=1)
            except:
                raise ValueError(f"Error concatenating {path}; check sample names. Sample names should be unique")
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
                lineage = ncbi.get_lineage(int(taxid))
                ranks = ncbi.get_rank(lineage)
                for tid, rank in ranks.items():
                    if rank == level:
                        taxid_map[taxid] = str(tid)
                        break
                if taxid not in taxid_map:  # If no match found at specified level
                    taxid_map[taxid] = taxid
            except (ValueError, KeyError):
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
            lineage = ncbi.get_lineage(int(j))
            ranks = ncbi.get_rank(lineage)
            for taxid, rank in ranks.items():
                if rank == i:
                    return str(taxid)
            return None
        except (ValueError, KeyError):
            return None

    def export(self, file: Optional[str] = None) -> pd.DataFrame:
        """Export annotations to CSV file.
        
        Args:
            file: Optional path to save CSV file
            
        Returns:
            DataFrame with exported annotations
        """
        df_return = self.DataFrame.loc[:, [col for col in self.DataFrame if col.startswith('taxID')]]
        lencol = [col for col in self.DataFrame if col.startswith('len')][0]
        df_return['length'] = self.DataFrame.loc[:, lencol].to_list()
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
        elif arg == "krakenunique":
            return "krakenunique"
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