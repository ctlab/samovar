"""Example ``samovar tools import --type table-scoring`` plugin.

KS test on pairwise Bray–Curtis distances. This is not the annotation
``--type scoring`` hook (``score(inputs, dest, config)``).
"""

from samovar.table_scorers import score_bray_ks_annotator


def score_annotator(annotation, tables_by_mode, annotator, modes, config=None):
    _ = config
    block = score_bray_ks_annotator(annotation, tables_by_mode, annotator, modes)
    block["scorer"] = "bray_ks_plugin"
    return block
