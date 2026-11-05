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
