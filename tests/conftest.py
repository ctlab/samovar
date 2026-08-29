from pathlib import Path

import pytest


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
        help="Import group: annotator, table_reads_generator, table_scoring, ...",
    )
