from pathlib import Path

import pytest

_CHECKOUT_POINTER = Path(__file__).resolve().parents[1] / "build" / "config_path"


@pytest.fixture(autouse=True)
def checkout_config_pointer_unchanged():
    """A test must not retarget the checkout install pointer."""
    before = _CHECKOUT_POINTER.read_text(encoding="utf-8") if _CHECKOUT_POINTER.is_file() else None
    yield
    after = _CHECKOUT_POINTER.read_text(encoding="utf-8") if _CHECKOUT_POINTER.is_file() else None
    if after == before:
        return
    if before is None:
        _CHECKOUT_POINTER.unlink(missing_ok=True)
    else:
        _CHECKOUT_POINTER.write_text(before, encoding="utf-8")
    raise AssertionError("test changed build/config_path")


def pytest_collection_modifyitems(config, items):
    """Unmarked tests are mandatory. CI runs ``pytest -m mandatory``.

    Mandatory is the core install path: Kraken2, Kaiju, and the prepare/exec
    wiring those two use. Kraken, KrakenUniq, and slow feature suites
    (apply, merge, visualization, MAG) are marked optional.
    """
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
    parser.addoption(
        "--db-path",
        action="store",
        default=None,
        help="Database directory for an annotator contract (else SAMOVAR_CONTRACT_DB)",
    )
