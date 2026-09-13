from pathlib import Path

import pytest


def pytest_collection_modifyitems(config, items):
    """Unmarked tests are mandatory; CI runs ``pytest -m mandatory``."""
    for item in items:
        names = {mark.name for mark in item.iter_markers()}
        if "optional" in names or "mandatory" in names:
            continue
        item.add_marker(pytest.mark.mandatory)


def pytest_addoption(parser):
    parser.addoption(
        "--tool",
        action="store",
        default=None,
        help="Path to a custom tool for tests/test_tool_contracts.py",
    )
    parser.addoption(
        "--tool-type",
        action="store",
        default=None,
        help="Import group: annotator, table_reads_generator, table_scoring, sample_scoring, export, ...",
    )
