"""Tests for the C++ sort-merge annotation combiner."""

import csv
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd
import pytest

from samovar.parse_annotators import Annotation, extract_true_taxid
from samovar.combine_tables import (
    apply_rank_level_csv,
    apply_species_level_csv,
    combine_with_cpp,
    ensure_combine_binary,
)


REPO = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def combiner():
    return ensure_combine_binary()


def _write_kaiju(path: Path, rows):
    with path.open("w", encoding="utf-8") as handle:
        for seq, taxid in rows:
            classified = "C" if str(taxid) not in {"0", ""} else "U"
            handle.write(f"{classified}\t{seq}\t{taxid}\n")


def _write_kraken2(path: Path, rows):
    with path.open("w", encoding="utf-8") as handle:
        for seq, taxid, length in rows:
            classified = "C" if str(taxid) not in {"0", ""} else "U"
            name = "unclassified" if str(taxid) in {"0", ""} else f"org (taxid {taxid})"
            handle.write(f"{classified}\t{seq}\t{name}\t{length}|{length}\t0:1\n")


def _write_custom(path: Path, rows):
    with path.open("w", encoding="utf-8") as handle:
        for seq, taxid in rows:
            handle.write(f"{seq}\t{taxid}\n")


def test_ensure_combine_binary(combiner):
    assert combiner.exists()
    assert os.access(combiner, os.X_OK)


def test_sort_merge_outer_join(tmp_path, combiner):
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    _write_kaiju(
        reports / "s1_kaiju.kaiju.out",
        [
            ("readB|taxid:562|", "562"),
            ("readA|taxid:9606|", "9606"),
            ("readA|taxid:9606|", "999"),
        ],
    )
    _write_kraken2(
        reports / "s1_kraken2.kraken2.out",
        [
            ("readC|taxid:4932|", "4932", "150"),
            ("readA|taxid:9606|", "9606", "126"),
        ],
    )
    subprocess.check_call(
        [
            str(combiner),
            "-i",
            str(reports),
            "-o",
            str(out),
            "-s",
            "1",
            "--chunk-rows",
            "2",
        ]
    )
    df = pd.read_csv(out / "s1.annotation.csv")
    df["seq"] = df["seq"].astype(str)
    by_seq = df.set_index("seq")
    assert set(by_seq.index) == {
        "readA|taxid:9606|",
        "readB|taxid:562|",
        "readC|taxid:4932|",
    }
    assert str(int(by_seq.loc["readA|taxid:9606|", "taxID_kaiju_0"])) == "9606"
    assert str(int(by_seq.loc["readA|taxid:9606|", "taxID_kraken2_1"])) == "9606"
    assert str(int(by_seq.loc["readB|taxid:562|", "taxID_kaiju_0"])) == "562"
    assert int(by_seq.loc["readB|taxid:562|", "taxID_kraken2_1"]) == 0
    assert int(by_seq.loc["readC|taxid:4932|", "taxID_kaiju_0"]) == 0
    assert str(int(by_seq.loc["readA|taxid:9606|", "true"])) == "9606"
    assert str(int(by_seq.loc["readC|taxid:4932|", "true"])) == "4932"
    assert int(by_seq.loc["readA|taxid:9606|", "length"]) == 126


def test_true_taxid_prefix_without_token(tmp_path, combiner):
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    _write_kaiju(
        reports / "1_kaiju.kaiju.out",
        [("Scer.fna-NC_001134.8_0_0", "0")],
    )
    subprocess.check_call(
        [str(combiner), "-i", str(reports), "-o", str(out), "-s", "1"]
    )
    df = pd.read_csv(out / "1.annotation.csv")
    assert str(int(df.loc[0, "true"])) == "4932"
    assert extract_true_taxid("Scer.fna-NC_001134.8_0_0") == "4932"


def test_truth_table_overrides_header_and_cami_ids(tmp_path, combiner):
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    _write_custom(
        reports / "sample_0_dummy.custom_dummy.out",
        [("S0R0", "9606"), ("S0R1", "9606"), ("readA|taxid:562|", "562")],
    )
    truth = tmp_path / "reads_mapping.tsv"
    truth.write_text(
        "#anonymous_read_id\tgenome_id\ttax_id\tread_id\n"
        "S0R0/1\tASV153.8\t818\tNZ_JBDQBI010000009.1-9644/1\n"
        "S0R1/1\tASV149.0\t762982\tNZ_GL883880.1-270717/1\n"
        "readA|taxid:562|\tignored\t9606\treadA\n",
        encoding="utf-8",
    )
    subprocess.check_call(
        [
            str(combiner),
            "-i",
            str(reports),
            "-o",
            str(out),
            "-s",
            "1",
            "--truth-table",
            str(truth),
        ]
    )
    df = pd.read_csv(out / "sample_0.annotation.csv")
    by_seq = df.set_index("seq")
    assert str(int(by_seq.loc["S0R0", "true"])) == "818"
    assert str(int(by_seq.loc["S0R1", "true"])) == "762982"
    assert str(int(by_seq.loc["readA|taxid:562|", "true"])) == "9606"


def test_truth_table_prefix_from_fasta_headers(tmp_path, combiner):
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    _write_custom(
        reports / "1_dummy.custom_dummy.out",
        [("NC_000913.3-0_0", "9606")],
    )
    truth = tmp_path / "gt.tsv"
    truth.write_text("seq\ttaxid\nNC_000913.3\t562\n", encoding="utf-8")
    subprocess.check_call(
        [
            str(combiner),
            "-i",
            str(reports),
            "-o",
            str(out),
            "--truth-table",
            str(truth),
        ]
    )
    df = pd.read_csv(out / "1.annotation.csv")
    assert str(int(df.loc[0, "true"])) == "562"


def _true_series(csv_path: Path) -> pd.Series:
    df = pd.read_csv(csv_path)
    return (
        df.assign(seq=df["seq"].astype(str))
        .set_index("seq")["true"]
        .map(lambda v: str(int(float(v))) if pd.notna(v) and str(v).strip() not in {"", "nan"} else "")
    )


def test_parse_genome_and_truth_table_agree_on_test_genomes(tmp_path, combiner):
    from samovar.ground_truth import append_from_genome_dir, iter_genome_truth_rows
    from samovar.seqio import is_fasta_name

    genomes = REPO / "data" / "test_genomes"
    truth = tmp_path / "gt.tsv"
    n = append_from_genome_dir(genomes / "meta", truth)
    n += append_from_genome_dir(genomes / "host", truth)
    assert n >= 3
    rows = []
    for folder in (genomes / "meta", genomes / "host"):
        for path in sorted(folder.iterdir()):
            if not path.is_file() or not is_fasta_name(path.name, protein=False):
                continue
            for token, tax in iter_genome_truth_rows(path):
                rows.append((f"{token}_0_0", "9606"))
    reports = tmp_path / "reports"
    reports.mkdir()
    _write_custom(reports / "1_dummy.custom_dummy.out", rows)
    parse_dir = tmp_path / "parse"
    table_dir = tmp_path / "table"
    subprocess.check_call(
        [str(combiner), "-i", str(reports), "-o", str(parse_dir), "-s", "1"]
    )
    subprocess.check_call(
        [
            str(combiner),
            "-i",
            str(reports),
            "-o",
            str(table_dir),
            "-s",
            "1",
            "--truth-table",
            str(truth),
        ]
    )
    parsed = _true_series(parse_dir / "1.annotation.csv")
    tabled = _true_series(table_dir / "1.annotation.csv")
    assert not (parsed == "").any()
    pd.testing.assert_series_equal(parsed, tabled, check_names=False)


def test_combine_with_cpp_accepts_gzipped_truth_table(tmp_path):
    import gzip as gzip_mod

    reports = tmp_path / "reports"
    reports.mkdir()
    _write_custom(reports / "s_dummy.custom_dummy.out", [("S0R0", "9606")])
    raw = tmp_path / "map.tsv"
    truth_gz = tmp_path / "map.tsv.gz"
    raw.write_text(
        "#anonymous_read_id\tgenome_id\ttax_id\tread_id\nS0R0/1\tg\t818\tr/1\n"
    )
    with raw.open("rb") as inn, gzip_mod.open(truth_gz, "wb") as out:
        out.write(inn.read())
    out_dir = tmp_path / "ann"
    combine_with_cpp(str(reports), str(out_dir), 1, truth_table=str(truth_gz))
    df = pd.read_csv(out_dir / "s.annotation.csv")
    assert str(int(df.loc[0, "true"])) == "818"


def test_sample_name_keeps_underscores(tmp_path, combiner):
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    _write_kaiju(reports / "sample_0_kaiju-test.kaiju.out", [("r0|taxid:562|", "562")])
    _write_kaiju(reports / "sample_2_kaiju-test.kaiju.out", [("r2|taxid:818|", "818")])
    _write_kraken2(
        reports / "1_full_kraken2-test.kraken2.out",
        [("r1|taxid:9606|", "9606", "150")],
    )
    subprocess.check_call(
        [str(combiner), "-i", str(reports), "-o", str(out), "-s", "1"]
    )
    names = sorted(p.name for p in out.glob("*.annotation.csv"))
    assert names == [
        "1_full.annotation.csv",
        "sample_0.annotation.csv",
        "sample_2.annotation.csv",
    ]
    assert "562" in (out / "sample_0.annotation.csv").read_text()
    assert "818" in (out / "sample_2.annotation.csv").read_text()


def test_large_tables_chunked_merge(tmp_path, combiner):
    n = 20000
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    kaiju = reports / "big_kaiju.kaiju.out"
    kraken = reports / "big_kraken2.kraken2.out"
    extra = "only_kraken|taxid:4932|"
    with kaiju.open("w", encoding="utf-8") as k, kraken.open("w", encoding="utf-8") as r:
        for i in range(n, 0, -1):
            seq = f"r{i:06d}|taxid:{562 if i % 2 == 0 else 9606}|"
            k.write(f"C\t{seq}\t{562 if i % 2 == 0 else 9606}\n")
            if i % 3 != 0:
                r.write(
                    f"C\t{seq}\torg (taxid {9606 if i % 2 else 562})\t100|100\t0:1\n"
                )
        r.write(f"C\t{extra}\torg (taxid 4932)\t80|80\t0:1\n")
    subprocess.check_call(
        [
            str(combiner),
            "-i",
            str(reports),
            "-o",
            str(out),
            "-s",
            "1",
            "--chunk-rows",
            "1023",
        ]
    )
    df = pd.read_csv(out / "big.annotation.csv")
    assert len(df) == n + 1
    assert df["seq"].is_monotonic_increasing
    assert df["seq"].is_unique
    only = df[df["seq"] == extra].iloc[0]
    assert int(only["taxID_kaiju_0"]) == 0
    assert int(only["taxID_kraken2_1"]) == 4932
    assert int(only["true"]) == 4932


def test_python_wrapper_species_pass(tmp_path):
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    _write_kaiju(reports / "1_kaiju.kaiju.out", [("seq|taxid:511145|", "511145")])
    combine_with_cpp(str(reports), str(out), split_n=1, chunk_rows=1000)
    csv_path = out / "1.annotation.csv"
    apply_species_level_csv(str(csv_path), level="species")
    df = pd.read_csv(csv_path)
    mapped = str(int(df.loc[0, "taxID_kaiju_0"]))
    assert mapped in {"511145", "562"}


def test_python_wrapper_genus_helper_does_not_change_combine_output(tmp_path):
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    _write_kaiju(reports / "1_kaiju.kaiju.out", [("seq|taxid:511145|", "511145")])
    combine_with_cpp(str(reports), str(out), split_n=1, chunk_rows=1000)
    csv_path = out / "1.annotation.csv"
    df = pd.read_csv(csv_path)
    assert str(int(df.loc[0, "taxID_kaiju_0"])) == "511145"
    assert str(int(df.loc[0, "true"])) == "511145"
    copy = tmp_path / "viz.csv"
    copy.write_text(csv_path.read_text())
    apply_rank_level_csv(str(copy), level="genus", cache_path=str(tmp_path / "taxid_genera_map.tsv"))
    mapped = pd.read_csv(copy)
    assert str(int(mapped.loc[0, "taxID_kaiju_0"])) == "561"
    assert str(int(pd.read_csv(csv_path).loc[0, "taxID_kaiju_0"])) == "511145"
    cache = (tmp_path / "taxid_genera_map.tsv").read_text()
    assert "taxid|genera_taxid" in cache
    assert "511145|561" in cache


def test_matches_pandas_annotation_on_fixture_logs(tmp_path, combiner):
    data = REPO / "tests" / "data"
    reports = tmp_path / "reports"
    out = tmp_path / "ann"
    reports.mkdir()
    shutil.copy(data / "kaiju.log", reports / "1_kaiju.kaiju.out")
    shutil.copy(data / "kraken2.log", reports / "1_kraken2.kraken2.out")
    subprocess.check_call(
        [str(combiner), "-i", str(reports), "-o", str(out), "-s", "1"]
    )
    cpp = pd.read_csv(out / "1.annotation.csv").set_index("seq")
    py = Annotation(
        {
            str(reports / "1_kaiju.kaiju.out"): "kaiju",
            str(reports / "1_kraken2.kraken2.out"): "kraken2",
        },
        "(?<=taxid:)[0-9]+",
    ).export()
    assert set(cpp.index) == set(py.index)
    shared = cpp.join(py, lsuffix="_cpp", rsuffix="_py", how="inner")
    for seq in list(shared.index)[:50]:
        for col_cpp, col_py in (
            ("taxID_kaiju_0_cpp", "taxID_kaiju_0_py"),
            ("taxID_kraken2_1_cpp", "taxID_kraken2_1_py"),
        ):
            a = str(shared.loc[seq, col_cpp]).replace(".0", "")
            b = str(shared.loc[seq, col_py]).replace(".0", "")
            if a in {"", "nan"}:
                a = "0"
            if b in {"", "nan"}:
                b = "0"
            assert a == b, (seq, col_cpp, a, b)


def test_tenbac_reports_merge(tmp_path, combiner):
    src = REPO / "samovar_10bac" / "initial_reports"
    if not src.exists():
        pytest.skip("samovar_10bac/initial_reports is not present")
    out = tmp_path / "ann"
    subprocess.check_call(
        [str(combiner), "-i", str(src), "-o", str(out), "-s", "1", "--chunk-rows", "4096"]
    )
    csvs = sorted(out.glob("*.annotation.csv"))
    assert csvs
    total = 0
    for path in csvs:
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle))
        assert rows
        assert "true" in rows[0]
        tax_cols = [c for c in rows[0] if c.startswith("taxID_")]
        assert len(tax_cols) >= 2
        total += len(rows)
        zeros = sum(1 for row in rows if (row.get("true") or "") == "0")
        assert zeros == 0
    assert total > 1000
