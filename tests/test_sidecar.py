from pathlib import Path

import pytest

from samovar import sidecar




@pytest.mark.parametrize("name", ["nanosim", "kraken2", "kaiju"])
def test_sidecar_health_checks_binary_modules_and_versions(tmp_path, monkeypatch, name):
    prefix = Path(tmp_path)
    binary = prefix / "bin" / sidecar.SIDECARS[name]["binary"]
    binary.parent.mkdir()
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)

    for helper in sidecar.SIDECARS[name].get("companions", ()):
        path = binary.parent / helper
        path.write_text("#!/bin/sh\n")
        path.chmod(0o755)
    monkeypatch.setattr(sidecar, "env_has_runtime", lambda *_: True)
    monkeypatch.setattr(sidecar, "env_has_versions", lambda *_: True)
    assert sidecar.sidecar_is_healthy(name, prefix)

    monkeypatch.setattr(sidecar, "env_has_versions", lambda *_: False)
    assert not sidecar.sidecar_is_healthy(name, prefix)


def test_motus_install_preserves_independently_configured_companions(tmp_path, monkeypatch):
    from samovar.main_config import lookup_tool_record
    existing = {'exec': {'env': '', 'path': str(tmp_path/'custom-samtools')}, 'type': 'aligner'}
    cfg = {'tools': {'samtools': existing.copy()}}
    monkeypatch.setattr(sidecar, 'load_config', lambda: cfg)
    monkeypatch.setattr(sidecar, 'write_config', lambda value: None)
    prefix = tmp_path/'motus-env'
    sidecar.record_sidecar('motus', prefix, str(prefix/'bin/motus'))
    assert cfg['tools']['samtools'] == existing
    assert lookup_tool_record(cfg, 'motus')['type'] == 'profiler'
    assert lookup_tool_record(cfg, 'bwa')['exec']['path'] == str(prefix/'bin/bwa')
    assert lookup_tool_record(cfg, 'bwa')['type'] == 'aligner'
