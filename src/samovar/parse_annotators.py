"""Module for parsing and processing taxonomic annotation outputs from various tools.

This module provides functionality to read and process taxonomic annotation outputs
from tools like Kraken, Kaiju, and MetaPhlAn. It standardizes the output format
to allow for easy comparison and analysis of results from different tools.
"""
import pandas as pd
import os
import re
from typing import Dict, List, Optional
from samovar.annotators_wrapper import get_annotator_instance
from samovar.taxonomy_engine import NCBITaxonomyParser

class Annotation:
    """Class for handling taxonomic annotations from multiple tools.
    
    Provides functionality to combine, process, and analyze results from different 
    classification tools at various taxonomic ranks.
    """
    def __init__(self, file_path_type: Dict[str, str], nodes_dmp_path: str, get_true_annotation: Optional[str] = None):
        """Initialize the Annotation object and combine input files.
        
        Args:
            file_path_type: Dictionary mapping file paths to their tool types (e.g., {'path/to/file': 'kraken2'}).
            nodes_dmp_path: Path to the NCBI nodes.dmp file for taxonomy resolution.
            get_true_annotation: Optional regex pattern to extract ground truth taxIDs from sequence headers.
        """
        self.id = 0
        self.DataFrame = pd.DataFrame()
        self.taxonomy = NCBITaxonomyParser(nodes_dmp_path)
        
        # Iteratively process each annotation file
        for path, tool_type in file_path_type.items():
            # Instantiate the appropriate annotator wrapper FIRST
            annotator = get_annotator_instance(tool_type, run_config={}, config={})
            
            # Concatenate the processed data into the main DataFrame
            try:
                # Check if file is empty before processing
                if os.path.getsize(path) == 0:
                    print(f"[WARNING] File is empty: {path}. Emulating unclassified predictions.")
                    df = pd.DataFrame(columns=["seq", "taxID"])
                else:
                    # Parse the raw output into a standard DataFrame [seq, taxID]
                    df = annotator.parse_output(path)
                
                # [FIX]: Do not skip empty DataFrames! We must concat them to create a column of 0s.
                if df.empty and os.path.getsize(path) != 0:
                    print(f"[WARNING] No data found in {path}. Treating entire sample as unclassified for this tool.")
            
                # Set read IDs as index and ensure taxIDs are strings for consistency
                df = df.set_index("seq").astype({"taxID": "string"})
                
                # Rename columns to include tool type and unique instance ID to avoid collisions
                df.columns = [f"{col}_{tool_type}_{self.id}" for col in df.columns]
                
                # Handle potential duplicate sequence IDs in the tool's output
                if df.index.duplicated().any():
                    df = df[~df.index.duplicated(keep="first")]
                
                # Concatenate the processed data into the main DataFrame
                try:
                    self.DataFrame = pd.concat([self.DataFrame, df], axis=1)
                except Exception as e:
                    raise ValueError(f"Failed to concatenate data from {path}: {e}")
                
                # Increment the ID for the next tool in the dictionary
                self.id += 1

            except Exception as e:
                print(f"[WARNING] Could not parse {path}: {e}. Skipping this file.")
                continue

        self.DataFrame = self.DataFrame.fillna("0")

        # [FIX]: Safer ground truth extraction using standard re.search to avoid pandas regex edge-cases
        self.true_annotation = []
        if get_true_annotation is not None:
            pattern = re.compile(get_true_annotation)
            for idx in self.DataFrame.index:
                match = pattern.search(str(idx))
                self.true_annotation.append(match.group(0) if match else "")
            print(f"Extracted {len([x for x in self.true_annotation if x])} ground truth annotations.")

        # Build unique lists of annotations and ranks for further analysis
        taxid_columns = self.tr()
        all_found_taxids = []
        for name, column in taxid_columns.items():
            all_found_taxids += list(set(column.dropna()))

        self.annotation_list = self.list2set(all_found_taxids)
        self.true_annotation_list = self.list2set(self.true_annotation)
        self.rank_list = self.list2set([*self.annotation_list, *self.true_annotation_list])

    def true_annotation_unique(self) -> set:
        """Return a set of unique ground truth annotations."""
        return set(self.true_annotation)

    def true_annotation_rename(self, change_dict: Dict[str, str]):
        """Rename ground truth labels using a mapping dictionary."""
        self.true_annotation = [change_dict.get(ta, "") for ta in self.true_annotation]

    def correct_level(self, level: str = "species") -> None:
        """Adjust all taxIDs in the DataFrame to the specified taxonomic level."""
        taxid_cols = [col for col in self.DataFrame.columns if col.startswith("taxID")]
        
        unique_taxids = set()
        for col in taxid_cols:
            unique_taxids.update(self.DataFrame[col].dropna().unique())
        
        taxid_map = {}
        for taxid in unique_taxids:
            if taxid == "0" or pd.isna(taxid):
                taxid_map[taxid] = taxid
                continue
            
            try:
                ancestor = self.taxonomy.get_ancestor_by_rank(int(taxid), level)
                taxid_map[taxid] = str(ancestor) if ancestor is not None else taxid
            except ValueError:
                taxid_map[taxid] = taxid
        
        for col in taxid_cols:
            self.DataFrame[col] = self.DataFrame[col].map(taxid_map)

    
    def get_taxid_at_rank(self, tax_id: str, target_rank: str) -> Optional[str]:
        """Retrieve the taxID at a specific rank for a given input taxID."""
        if tax_id == "0" or pd.isna(tax_id):
            return "0"
        try:
            ancestor = self.taxonomy.get_ancestor_by_rank(int(tax_id), target_rank)
            return str(ancestor) if ancestor is not None else None
        except ValueError:
            return None

    def list2rank(self, taxid_list: List, target_rank: str) -> List[Optional[str]]:
        """Convert a list of taxIDs to their corresponding IDs at a target rank."""
        unique_taxids = self.list2set(taxid_list)
        return [self.get_taxid_at_rank(i, target_rank) for i in unique_taxids]

    def rank_annotation(self, rank: str = "species") -> "RankAnnotation":
        """Map all current annotations to a specific taxonomic rank."""
        rank_list = [self.get_taxid_at_rank(j, rank) for j in self.rank_list]
        rank_dict = dict(zip(self.rank_list, rank_list))
        return RankAnnotation(rank).make(self.full(), rank_dict)

    def expand_annotation(self, rank: List[str] = ["species"]) -> "ExpandAnnotation":
        """Generate annotations for multiple taxonomic ranks simultaneously."""
        full_rank_annotation = ExpandAnnotation()
        for r in rank:
            full_rank_annotation.add(self.rank_annotation(r))
        return full_rank_annotation

    def correct_annotations(self, rank: str = "species") -> pd.DataFrame:
        """Calculate counts of correct predictions at a specific rank."""
        return pd.DataFrame(self.rank_annotation(rank).correct_annotation().annotation.value_counts())

    def full(self) -> pd.DataFrame:
        """Return the complete DataFrame including the ground truth column."""
        tmp = self.tr()
        tmp["true"] = self.true_annotation
        return tmp

    def tr(self) -> pd.DataFrame:
        """Filter the DataFrame to include only taxID columns."""
        return self.DataFrame.copy().filter(regex="taxID.*")

    @staticmethod
    def list2set(a: List) -> List[str]:
        """Convert a list to a unique set of strings."""
        return list(set([str(i) for i in a]))

    def export(self, file: Optional[str] = None) -> pd.DataFrame:
        """Export the taxID columns and ground truth to a TSV file."""
        df_return = self.DataFrame.loc[:, [col for col in self.DataFrame if col.startswith("taxID")]]
        df_return["true"] = self.true_annotation
        
        if file is not None:
            df_return.to_csv(file, sep="\t", index=True) 
        return df_return


class RankAnnotation:
    """Class representing taxonomic annotations at a single fixed rank."""
    
    def __init__(self, rank: str):
        self.rank = rank
        self.annotation = pd.DataFrame()

    def add(self, name: str, annotation: List):
        self.annotation[name] = annotation

    def make(self, annotation: pd.DataFrame, rank_dict: Dict[str, str]) -> "RankAnnotation":
        for name, column in annotation.items():
            self.add(str(name), [rank_dict.get(str(ta)) for ta in column])
        self.reindex(annotation.index)
        return self

    def reindex(self, index):
        self.annotation.index = index

    def y(self) -> pd.Series:
        return self.annotation["true"]

    def x(self) -> pd.DataFrame:
        return self.annotation.copy().drop("true", axis=1)

    def correct_annotation(self) -> "RankAnnotation":
        tmp = RankAnnotation(self.rank)
        for name, column in self.x().items():
            tmp.add(str(name), pd.DataFrame(column == self.y()))
        tmp.reindex(self.annotation.index)
        return tmp


class ExpandAnnotation:
    """Class managing taxonomic annotations across multiple ranks."""
    
    def __init__(self):
        self.rank_annotation = {}

    def add(self, rank_annotation: RankAnnotation):
        self.rank_annotation[rank_annotation.rank] = rank_annotation

    def get(self, rank: str) -> RankAnnotation:
        return self.rank_annotation[rank]


def match_annotation(annotation_name: str) -> Optional[str]:
    """Determine the tool type from a filename (expects .tool_name.out format)."""
    basename = os.path.basename(annotation_name)
    if basename.endswith(".out"):
        parts = basename.split(".")
        if len(parts) >= 2:
            tool = parts[-2].lower()
            if tool.startswith("custom_"):
                return tool[7:]
            return tool
    return None

if __name__ == "__main__":
    import glob
    import sys

    # --- Configuration ---
    results_dir = "results"
    nodes_dmp = "db/kraken2/taxonomy/nodes.dmp"
    output_path = "results/merged_table.tsv"
    truth_regex = r"\|(\d+)" 

    print(f"[INFO] Searching for annotation files in {results_dir}...")
    out_files = glob.glob(os.path.join(results_dir, "*.out"))
    
    if not out_files:
        print(f"[ERROR] No .out files found in {results_dir}. Did Snakemake finish successfully?")
        sys.exit(1)

    file_mapping = {}
    for f in out_files:
        try:
            tool_type = match_annotation(f)
            if tool_type:
                file_mapping[f] = tool_type
                print(f"[INFO] Found {tool_type} output: {os.path.basename(f)}")
        except ValueError as e:
            print(f"[WARNING] Skipping file: {e}")

    if not file_mapping:
        print("[ERROR] Could not identify any tool types from the filenames.")
        sys.exit(1)

    print(f"[INFO] Initializing Annotation engine with {len(file_mapping)} files...")
    try:
        ann = Annotation(
            file_path_type=file_mapping, 
            nodes_dmp_path=nodes_dmp,
            get_true_annotation=truth_regex 
        )
        
        print(f"[INFO] Exporting combined data to {output_path}...")
        ann.export(output_path)
        print(f"[SUCCESS] Merged table created successfully at {output_path}")

    except Exception as e:
        print(f"[ERROR] An error occurred during processing: {e}")
