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

    Mandatory covers pipeline behavior and extension contracts with small local
    fixtures and core tools. Optional covers native ecosystem integrations,
    network downloads, detailed rendering and resource/algorithm tuning.
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


@pytest.fixture(scope="session")
def core_runtime():
    """Real core executable paths, without a developer's custom tool catalog."""
    import sys
    from samovar.paths import load_config
    from samovar.main_config import iter_tool_records

    installed = load_config()
    core = {name: row for name, row in iter_tool_records(installed).items()
            if name.split(":", 1)[0].lower() in {"iss", "snakemake", "kraken2", "kaiju"}}
    return {"root": str(_CHECKOUT_POINTER.parents[1]),
            "compilers": {**installed.get("compilers", {}), "python": sys.executable},
            "tools": core, "API": {"ncbi_email": "test@samovar.com"}}


@pytest.fixture(autouse=True)
def mandatory_environment(request, tmp_path, monkeypatch, core_runtime):
    """Core tests use private config/cache and a tiny taxonomy, never live NCBI."""
    if request.node.get_closest_marker("optional") is not None:
        return
    import json
    import socket
    import urllib.request
    from urllib.parse import urlsplit

    dump = Path(__file__).parent / "data" / "taxdump"
    config = {**core_runtime, "genomes": {"taxdump": str(dump.resolve())}}
    target = tmp_path / "install.json"
    target.write_text(json.dumps(config))
    monkeypatch.setenv("SAMOVAR_CONFIG", str(target))
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    monkeypatch.setenv("SAMOVAR_DATABASE", str(tmp_path / "database"))
    monkeypatch.setenv("SAMOVAR_PROCESSED_GENOMES", str(tmp_path / "processed"))
    monkeypatch.setenv("SAMOVAR_TAXDUMP", str(dump.resolve()))
    monkeypatch.setenv("SAMOVAR_OFFLINE", "1")

    def local_only(original):
        def checked(address, *args, **kwargs):
            host = address[0] if isinstance(address, tuple) else address
            if str(host) not in {"127.0.0.1", "::1", "localhost"}:
                pytest.fail(f"mandatory test attempted external network: {host}; "
                            "use a local fixture or mark the real integration optional")
            return original(address, *args, **kwargs)
        return checked

    monkeypatch.setattr(socket, "getaddrinfo", local_only(socket.getaddrinfo))
    monkeypatch.setattr(socket, "create_connection", local_only(socket.create_connection))
    # A local HTTP proxy can otherwise tunnel external traffic past socket checks.
    real_open = urllib.request.OpenerDirector.open

    def open_local(opener, request, *args, **kwargs):
        url = getattr(request, "full_url", request)
        host = urlsplit(str(url)).hostname
        if host and host not in {"127.0.0.1", "::1", "localhost"}:
            pytest.fail(f"mandatory test attempted external URL: {url}; "
                        "mock the transport or run the integration as optional")
        return real_open(opener, request, *args, **kwargs)

    monkeypatch.setattr(urllib.request.OpenerDirector, "open", open_local)
