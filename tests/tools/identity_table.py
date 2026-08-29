"""Identity table regenerator for contract tests (``--type table``)."""

from samovar.abundance import input_to_abundance_tables, normalize_abundance_table
from samovar.regenerate import cap_abundance_table, max_genomes_from_config


def regenerate(data, metadata, config):
    _ = metadata
    limit = max_genomes_from_config(config)
    return {
        name: cap_abundance_table(normalize_abundance_table(frame), limit)
        for name, frame in input_to_abundance_tables(data).items()
    }
