import os
from typing import Dict, Optional, Tuple


class NCBITaxonomyParser:
    """
    A lightweight, memory-efficient parser for the NCBI nodes.dmp file.
    Replaces the heavy ete3 dependency by loading only necessary
    parent-child relationships and taxonomic ranks into a hash map.
    """

    def __init__(self, nodes_dmp_path: str):
        # Dictionary to store the tree.
        # Key: tax_id (int)
        # Value: tuple (parent_tax_id: int, rank: str)
        self._tree: Dict[int, Tuple[int, str]] = {}
        self._load_tree(nodes_dmp_path)

    def _load_tree(self, path: str) -> None:
        """Parses the nodes.dmp file and populates the internal dictionary."""
        if not os.path.exists(path):
            raise FileNotFoundError(f"Taxonomy file not found: {path}")

        # Open the file and read line by line to prevent RAM overload
        with open(path, "r") as file:
            for line in file:
                # NCBI dump files use "\t|\t" as a column separator
                parts = line.split("\t|\t")

                # Ensure the line has enough columns
                if len(parts) >= 3:
                    tax_id = int(parts[0].strip())
                    parent_id = int(parts[1].strip())
                    rank = parts[2].strip()

                    # Store as a lightweight tuple
                    self._tree[tax_id] = (parent_id, rank)

    def get_lineage_ranks(self, tax_id: int) -> Dict[str, int]:
        """
        Traverses the taxonomy tree upwards from the given tax_id to the root.
        Returns a dictionary mapping taxonomic ranks to their respective TaxIDs.

        Example:
        get_lineage_ranks(562) -> {'species': 562, 'genus': 561, 'family': 543}
        """
        lineage = {}
        current_id = tax_id

        # Traverse until we hit the root (TaxID 1) or an unknown ID
        while current_id in self._tree and current_id != 1:
            parent_id, rank = self._tree[current_id]

            # Save the rank and ID (ignore uninformative "no rank" nodes)
            if rank != "no rank":
                lineage[rank] = current_id

            # Break circular references if the database is corrupted
            if parent_id == current_id:
                break

            # Move exactly one step up the phylogenetic tree
            current_id = parent_id

        return lineage

    def get_ancestor_by_rank(self, tax_id: int, target_rank: str) -> Optional[int]:
        """
        Returns the ancestor TaxID for a specific taxonomic rank, or None.
        """
        lineage = self.get_lineage_ranks(tax_id)
        return lineage.get(target_rank)
