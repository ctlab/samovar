import os
import subprocess
import tempfile
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
import yaml

from samovar.table2iss import (
    _resolve_r_executable,
    process_abundance_table,
    samovar_annotation_regenerate,
)


def _r_with_samovar_available() -> bool:
    """True when R can load the samovaR package (real integration path)."""
    try:
        r_path, r_lib_path = _resolve_r_executable()
    except FileNotFoundError:
        return False
    env = os.environ.copy()
    if r_lib_path:
        env["R_LIBS"] = r_lib_path
        env["R_LIBS_USER"] = r_lib_path
    try:
        proc = subprocess.run(
            [r_path, "--vanilla", "-s", "-e", "library(samovaR)"],
            capture_output=True,
            text=True,
            timeout=120,
            env=env,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    return proc.returncode == 0


requires_r_samovar = pytest.mark.skipif(
    not _r_with_samovar_available(),
    reason="R with samovaR package is required for this integration test",
)


@pytest.fixture
def test_data_dir():
    return Path("data/test_annotations")


@pytest.fixture
def test_output_dir():
    output_dir = Path("tests_outs/test_samovar")
    output_dir.mkdir(parents=True, exist_ok=True)
    return output_dir


@pytest.fixture
def mock_config(test_output_dir):
    return {
        "threshold_amount": 1e-5,
        "plot_log": False,
        "min_cluster_size": 2,
        "N": 10,
        "N_reads": 1000,
        "output_dir": str(test_output_dir),
    }


def _write_config(mock_config):
    handle = tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False)
    yaml.dump(mock_config, handle)
    handle.close()
    return handle.name


def test_resolve_r_executable_uses_path_when_config_missing(tmp_path, monkeypatch):
    monkeypatch.setattr(
        "samovar.table2iss.os.path.dirname",
        lambda _path: str(tmp_path),
    )
    monkeypatch.setattr(
        "samovar.table2iss.shutil.which",
        lambda name: "/usr/bin/R" if name == "R" else None,
    )
    resolved, lib = _resolve_r_executable()
    assert resolved == "/usr/bin/R"
    assert lib is None


def test_resolve_r_executable_raises_when_missing(monkeypatch, tmp_path):
    monkeypatch.setattr(
        "samovar.table2iss.os.path.dirname",
        lambda _path: str(tmp_path),
    )
    monkeypatch.setattr("samovar.table2iss.shutil.which", lambda _name: None)
    with pytest.raises(FileNotFoundError, match="R executable not found"):
        _resolve_r_executable()


def test_samovar_annotation_regenerate_invokes_r(tmp_path, mock_config):
    """Unit test: wrapper builds the R command and fails if R exits non-zero."""
    annotation_dir = tmp_path / "ann"
    annotation_dir.mkdir()
    output_dir = tmp_path / "out"
    config_path = _write_config({**mock_config, "output_dir": str(output_dir)})

    calls = []

    def fake_run(cmd, check=True, env=None):
        calls.append(cmd)
        class Result:
            returncode = 0

        # Simulate R writing an abundance CSV
        output_dir.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"taxid": [562], "N_1": [10]}).to_csv(
            output_dir / "kaiju.csv", index=False
        )
        return Result()

    try:
        with patch("samovar.table2iss._resolve_r_executable", return_value=("R", None)), \
             patch("samovar.table2iss.subprocess.run", side_effect=fake_run):
            samovar_annotation_regenerate(
                annotation_dir=str(annotation_dir),
                config_samovar=config_path,
                output_dir=str(output_dir),
            )
    finally:
        os.unlink(config_path)

    assert len(calls) == 1
    cmd = calls[0]
    assert cmd[0] == "R"
    assert "--config" in cmd
    assert "--annotation_dir" in cmd
    assert "--output_dir" in cmd
    assert str(annotation_dir) in cmd
    assert str(output_dir) in cmd
    assert (output_dir / "kaiju.csv").exists()


def test_samovar_annotation_regenerate_raises_on_r_failure(tmp_path, mock_config):
    annotation_dir = tmp_path / "ann"
    annotation_dir.mkdir()
    output_dir = tmp_path / "out"
    config_path = _write_config({**mock_config, "output_dir": str(output_dir)})

    def boom(cmd, check=True, env=None):
        raise subprocess.CalledProcessError(returncode=1, cmd=cmd)

    try:
        with patch("samovar.table2iss._resolve_r_executable", return_value=("R", None)), \
             patch("samovar.table2iss.subprocess.run", side_effect=boom):
            with pytest.raises(RuntimeError, match="failed while running R"):
                samovar_annotation_regenerate(
                    annotation_dir=str(annotation_dir),
                    config_samovar=config_path,
                    output_dir=str(output_dir),
                )
    finally:
        os.unlink(config_path)


@requires_r_samovar
def test_samovar_annotation_regenerate_basic(test_data_dir, test_output_dir, mock_config):
    """Integration: real R + samovaR regeneration produces abundance CSVs."""
    assert test_data_dir.exists(), f"Missing fixture annotations at {test_data_dir}"
    config_path = _write_config(mock_config)

    try:
        for stale in test_output_dir.glob("*.csv"):
            stale.unlink()

        samovar_annotation_regenerate(
            annotation_dir=str(test_data_dir),
            config_samovar=config_path,
            output_dir=str(test_output_dir),
        )

        output_files = list(test_output_dir.glob("*.csv"))
        assert len(output_files) > 0, "No output files were created"

        for output_file in output_files:
            df = pd.read_csv(output_file)
            assert not df.empty, f"Output file {output_file} is empty"
            assert "taxid" in df.columns
            numeric = df.drop(columns=["taxid"], errors="ignore")
            if not numeric.empty:
                values = numeric.apply(pd.to_numeric, errors="coerce")
                assert (values.fillna(0) >= 0).all().all()
    finally:
        os.unlink(config_path)


@requires_r_samovar
def test_samovar_annotation_regenerate_integration(test_data_dir, test_output_dir, mock_config):
    """Integration: regenerated tables feed process_abundance_table."""
    config_path = _write_config(mock_config)

    try:
        for stale in test_output_dir.glob("*.csv"):
            stale.unlink()

        samovar_annotation_regenerate(
            annotation_dir=str(test_data_dir),
            config_samovar=config_path,
            output_dir=str(test_output_dir),
        )

        output_files = list(test_output_dir.glob("*.csv"))
        assert len(output_files) > 0, "No output files were created"
        test_file = output_files[0]

        with tempfile.TemporaryDirectory() as genome_dir:
            with patch("samovar.table2iss.fetch_genome") as mock_fetch, \
                 patch("samovar.table2iss.get_genome_file") as mock_get, \
                 patch("samovar.table2iss.regenerate_metagenome") as mock_regenerate:
                mock_get.return_value = None
                mock_fetch.return_value = os.path.join(genome_dir, "dummy.fasta")
                with open(os.path.join(genome_dir, "dummy.fasta"), "w") as f:
                    f.write(">dummy\nATCG\n")
                mock_regenerate.return_value = None

                result = process_abundance_table(
                    table=str(test_file),
                    genome_dir=genome_dir,
                    output_dir=str(test_output_dir / "processed"),
                )

                assert isinstance(result, pd.DataFrame)
                assert not result.empty
    finally:
        os.unlink(config_path)


def test_process_abundance_table_after_mocked_regenerate(tmp_path, mock_config):
    """CI-safe path: abundance CSV → process_abundance_table without R."""
    table = tmp_path / "kaiju.csv"
    pd.DataFrame({"taxid": [562, 9606], "N_1": [5, 3]}).to_csv(table, index=False)
    genome_dir = tmp_path / "genomes"
    genome_dir.mkdir()
    for taxid in (562, 9606):
        (genome_dir / f"{taxid}.fasta").write_text(f">{taxid}\nATCG\n")
    out = tmp_path / "processed"

    with patch("samovar.table2iss.regenerate_metagenome") as mock_regenerate:
        mock_regenerate.return_value = None
        result = process_abundance_table(
            table=str(table),
            genome_dir=str(genome_dir),
            output_dir=str(out),
        )

    assert isinstance(result, pd.DataFrame)
    assert not result.empty
    assert set(result["taxid"].astype(str)) == {"562", "9606"}
    assert mock_regenerate.called
