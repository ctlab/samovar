"""Identity table regenerator for contract tests (``--type table``)."""

from samovar.abundance import input_to_abundance_tables, normalize_abundance_table


def regenerate(data, metadata, config):
    _ = metadata, config
    return {
        name: normalize_abundance_table(frame)
        for name, frame in input_to_abundance_tables(data).items()
    }
