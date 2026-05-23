import os
from typing import Dict, Optional, Tuple


class NCBITaxonomyParser:
    """Lightweight parser for NCBI ``nodes.dmp``.

    Stores parent/rank as a hash map instead of loading a full ete3 database.
    Used as a fallback (and as a standalone helper) when resolving taxIDs to
    a rank such as genus or species.
    """

    def __init__(self, nodes_dmp_path: str):
        self.nodes_dmp_path = nodes_dmp_path
        self._tree: Dict[int, Tuple[int, str]] = {}
        self._load_tree(nodes_dmp_path)

    @property
    def tree(self) -> Dict[int, Tuple[int, str]]:
        return self._tree

    def _load_tree(self, path: str) -> None:
        if not os.path.exists(path):
            raise FileNotFoundError(f"Taxonomy file not found: {path}")

        with open(path, "r", encoding="utf-8", errors="ignore") as file:
            for line in file:
                parts = line.split("\t|\t")
                if len(parts) < 3:
                    continue
                try:
                    tax_id = int(parts[0].strip())
                    parent_id = int(parts[1].strip())
                except ValueError:
                    continue
                rank = parts[2].strip().rstrip("|").strip()
                self._tree[tax_id] = (parent_id, rank)

    def get_lineage_ranks(self, tax_id: int) -> Dict[str, int]:
        """Walk from ``tax_id`` to the root.

        Example: ``get_lineage_ranks(562) -> {'species': 562, 'genus': 561, ...}``
        """
        lineage: Dict[str, int] = {}
        current_id = tax_id

        while current_id in self._tree and current_id != 1:
            parent_id, rank = self._tree[current_id]
            if rank and rank != "no rank":
                lineage[rank] = current_id
            if parent_id == current_id:
                break
            current_id = parent_id

        return lineage

    def get_ancestor_by_rank(self, tax_id: int, target_rank: str) -> Optional[int]:
        """Return the ancestor (or self) at ``target_rank``, or None."""
        lineage = self.get_lineage_ranks(tax_id)
        return lineage.get(target_rank)
