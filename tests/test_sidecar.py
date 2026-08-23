from pathlib import Path

from samovar import sidecar


def test_nanosim_sidecar_pins_profile_compatible_runtime():
    spec = sidecar.SIDECARS["nanosim"]
    assert spec["python"] == "3.8"
    assert "numpy=1.23.5" in spec["packages"]
    assert "scikit-learn=0.23.2" in spec["packages"]
    assert "htseq" in spec["packages"]


def test_sidecar_health_checks_binary_modules_and_versions(tmp_path, monkeypatch):
    prefix = Path(tmp_path)
    binary = prefix / "bin" / "simulator.py"
    binary.parent.mkdir()
    binary.write_text("#!/bin/sh\n")
    binary.chmod(0o755)

    monkeypatch.setattr(sidecar, "env_has_runtime", lambda *_: True)
    monkeypatch.setattr(sidecar, "env_has_versions", lambda *_: True)
    assert sidecar.sidecar_is_healthy("nanosim", prefix)

    monkeypatch.setattr(sidecar, "env_has_versions", lambda *_: False)
    assert not sidecar.sidecar_is_healthy("nanosim", prefix)
