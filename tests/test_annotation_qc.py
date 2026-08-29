"""Empty annotator output: warn, skip table regen, fail when nothing remains."""

import pandas as pd
import pytest

from samovar.abundance import materialize_observed_abundance
from samovar.annotation_qc import (
    EmptyAnnotatorsError,
    diagnose_annotator,
    filter_classified_abundance_tables,
)
from samovar.annotators_wrapper import skip_empty_reads_cmd


def test_filter_drops_unclassified_keeps_classified(tmp_path):
    tables = {
        "kaiju": pd.DataFrame({"taxid": ["0"], "N_1": [10]}),
        "kraken2": pd.DataFrame({"taxid": ["562", "0"], "N_1": [4, 1]}),
    }
    with pytest.warns(UserWarning, match="kaiju"):
        kept = filter_classified_abundance_tables(tables, fatal_if_none=True)
    assert "kraken2" in kept
    assert "kaiju" not in kept


def test_filter_all_unclassified_is_error():
    tables = {
        "kaiju": pd.DataFrame({"taxid": ["0"], "N_1": [3]}),
        "kraken2": pd.DataFrame({"taxid": ["unclassified"], "N_1": [8]}),
    }
    with pytest.raises(EmptyAnnotatorsError, match="All annotators"):
        with pytest.warns(UserWarning):
            filter_classified_abundance_tables(tables, fatal_if_none=True)


def test_diagnose_did_not_run_vs_unclassified(tmp_path):
    reports = tmp_path / "initial_reports"
    reports.mkdir()
    (reports / "1_kaiju.kaiju.out").write_text("")
    (reports / "1_kraken2.kraken2.out").write_text("U\tread1\t0\n")
    empty = pd.DataFrame({"taxid": ["0"], "N_1": [1]})
    assert diagnose_annotator("kaiju", empty, reports_dir=reports) == "did_not_run"
    assert diagnose_annotator("kraken2", empty, reports_dir=reports) == "unclassified"
    msg_run = diagnose_annotator("kaiju", empty, reports_dir=reports)
    from samovar.annotation_qc import warning_message

    assert "did not write output" in warning_message("kaiju", msg_run)
    assert "no taxa in this community" in warning_message("kraken2", "unclassified")


def test_materialize_skips_empty_annotator(tmp_path):
    ann = tmp_path / "initial_annotations"
    ann.mkdir()
    pd.DataFrame(
        {
            "seq": ["a", "b"],
            "taxID_kaiju_0": [0, 0],
            "taxID_kraken2_0": [562, 562],
        }
    ).to_csv(ann / "1.annotation.csv", index=False)
    reports = tmp_path / "initial_reports"
    reports.mkdir()
    (reports / "1_kaiju.kaiju.out").write_text("")
    (reports / "1_kraken2.kraken2.out").write_text(
        "C\ta\torg (taxid 562)\t10|10\t0:1\n"
    )
    with pytest.warns(UserWarning, match="kaiju"):
        tables = materialize_observed_abundance(tmp_path)
    assert "kraken2" in tables
    assert "kaiju" not in tables
    written = {p.stem for p in (tmp_path / "initial_abundance").glob("*.csv")}
    assert "kraken2" in written
    assert "kaiju" not in written


def test_materialize_all_unclassified_raises(tmp_path):
    ann = tmp_path / "initial_annotations"
    ann.mkdir()
    pd.DataFrame(
        {"seq": ["a"], "taxID_kaiju_0": [0], "taxID_kraken2_0": [0]}
    ).to_csv(ann / "1.annotation.csv", index=False)
    with pytest.raises(EmptyAnnotatorsError):
        with pytest.warns(UserWarning):
            materialize_observed_abundance(tmp_path)


def test_skip_empty_reads_one_byte_stub(tmp_path):
    import os

    stub = tmp_path / "s_R1.fastq"
    stub.write_bytes(b"x")
    out = tmp_path / "s.out"
    cmd = skip_empty_reads_cmd(str(stub), [str(out)], "false")
    rc = os.system(cmd)
    assert rc == 0
    assert out.exists()
    assert "fastq_has_reads" in cmd
