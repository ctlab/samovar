from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pandas as pd

from samovar.cami_profile import (
    CAMI_FORMAT_VERSION,
    CAMI_RANKS,
    counts_to_cami_profile,
    dump_cami,
    write_cami_profile,
)
from samovar.kreport import KReportTaxonomy
from samovar.parse_annotators import Annotation
from samovar.taxonomy import GTDB_CAMI_RANKS
from tests.test_kreport import write_mini_taxdump

_HEADER_KV = re.compile(
    r"^@(_[A-Za-z]*_)?[A-Za-z]+[A-Za-z0-9]*:[A-Za-z0-9,.;_|]*$"
)
_PCT = re.compile(r"^[0-9]+(\.[0-9]{0,6})?$")
_SAMPLEID = re.compile(r"^[A-Za-z0-9._]+$")


def _sections(text: str):
    chunks = re.split(r"\n\s*\n", text.strip())
    return [c.strip() for c in chunks if c.strip()]


def parse_profile(text: str):
    sections = []
    for chunk in _sections(text):
        header = {}
        columns = None
        rows = []
        for line in chunk.splitlines():
            if not line.strip() or line.startswith("#"):
                continue
            if line.startswith("@@"):
                columns = line[2:].split("\t")
                continue
            if line.startswith("@"):
                assert _HEADER_KV.match(line), line
                key, value = line[1:].split(":", 1)
                header[key.lower()] = value
                continue
            fields = line.split("\t")
            assert columns is not None
            assert len(fields) == len(columns), (columns, fields)
            rows.append(dict(zip(columns, fields)))
        sections.append((header, columns, rows))
    return sections


def assert_rfc(text: str):
    sections = parse_profile(text)
    assert sections
    sample_ids = []
    for header, columns, rows in sections:
        assert header["sampleid"]
        assert _SAMPLEID.match(header["sampleid"])
        sample_ids.append(header["sampleid"])
        assert header["version"] == CAMI_FORMAT_VERSION
        tax_sys = header.get("taxonomyid", "ncbi").lower()
        assert tax_sys in {"ncbi", "gtdb"}
        expected_ranks = GTDB_CAMI_RANKS if tax_sys == "gtdb" else CAMI_RANKS
        assert header["ranks"] == "|".join(expected_ranks)
        assert columns[:4] == ["TAXID", "RANK", "TAXPATH", "PERCENTAGE"] or columns[
            :5
        ] == ["TAXID", "RANK", "TAXPATH", "TAXPATHSN", "PERCENTAGE"]
        assert columns[0] == "TAXID"
        assert columns[1] == "RANK"
        assert columns[2] == "TAXPATH"
        assert columns[-1] == "PERCENTAGE"
        rank_sums = defaultdict(float)
        by_id = {}
        for row in rows:
            taxid = row["TAXID"]
            rank = row["RANK"]
            path = row["TAXPATH"]
            pct = row["PERCENTAGE"]
            assert taxid and taxid != "0"
            assert rank in expected_ranks or rank == ""
            assert _PCT.match(pct), pct
            assert path.split("|")[-1] == str(taxid)
            if "TAXPATHSN" in row:
                assert len(row["TAXPATH"].split("|")) == len(row["TAXPATHSN"].split("|"))
            value = float(pct)
            assert 0 <= value <= 100
            if rank:
                rank_sums[rank] += value
            by_id[taxid] = row
        for rank, total in rank_sums.items():
            assert total <= 100.000001, (rank, total)
        for row in rows:
            parts = [p for p in row["TAXPATH"].split("|") if p]
            if len(parts) < 2:
                continue
            parent = parts[-2]
            if parent in by_id:
                assert float(by_id[parent]["PERCENTAGE"]) + 1e-9 >= float(row["PERCENTAGE"])
    assert len(sample_ids) == len(set(sample_ids))
    return sections


def test_cami_rfc_taxpath_and_percentages(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    tax = KReportTaxonomy.from_taxdump(write_mini_taxdump(tmp_path / "taxdump"))
    text = counts_to_cami_profile({562: 2, 0: 1}, tax, "s1", total=3)
    sections = assert_rfc(text)
    header, columns, rows = sections[0]
    assert header["sampleid"] == "s1"
    assert "TAXPATHSN" in columns
    by_id = {row["TAXID"]: row for row in rows}
    assert by_id["2"]["RANK"] == "superkingdom"
    assert by_id["2"]["TAXPATH"] == "2"
    assert by_id["2"]["TAXPATHSN"] == "Bacteria"
    assert by_id["561"]["RANK"] == "genus"
    assert by_id["561"]["TAXPATH"] == "2|||||561"
    assert by_id["562"]["RANK"] == "species"
    assert by_id["562"]["TAXPATH"] == "2|||||561|562"
    assert by_id["562"]["TAXPATHSN"].endswith("Escherichia coli")
    assert abs(float(by_id["562"]["PERCENTAGE"]) - 66.666667) < 1e-6
    assert "0" not in by_id
    strain = counts_to_cami_profile({83333: 2}, tax, "s2")
    strain_rows = parse_profile(strain)[0][2]
    leaf = next(r for r in strain_rows if r["TAXID"] == "83333")
    assert leaf["RANK"] == ""
    assert leaf["TAXPATH"].endswith("|83333")


def test_annotation_to_cami_and_dump(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    obj = Annotation.from_long_table(
        pd.DataFrame(
            {
                "seq": ["a", "b", "c"],
                "taxID_kaiju_0": ["562", "562", "0"],
                "sample": ["1", "1", "1"],
            }
        )
    )
    reports = obj.to_cami(taxdump=dump)
    assert_rfc(reports["kaiju"]["1"])
    dest = tmp_path / "out.profile"
    written = obj.to_CAMI(dest, taxdump=dump)
    assert Path(written).suffix == ".profile"
    assert_rfc(Path(written).read_text())


def test_dump_cami_concatenates_samples(tmp_path, monkeypatch):
    monkeypatch.setenv("SAMOVAR_CACHE_DIR", str(tmp_path / "cache"))
    dump = write_mini_taxdump(tmp_path / "taxdump")
    tables = {
        "kaiju": pd.DataFrame({"taxid": [562, 0], "N_1": [2, 1], "N_2": [1, 0]}),
    }
    dest = dump_cami(tables, tmp_path / "profiles", taxdump=dump)
    text = (dest / "kaiju.profile").read_text()
    sections = assert_rfc(text)
    assert {s[0]["sampleid"] for s in sections} == {"1", "2"}


def test_write_cami_profile_fallback_without_taxdump(tmp_path, monkeypatch):
    monkeypatch.setattr("samovar.cami_profile.try_taxonomy", lambda taxdump=None: None)
    path = write_cami_profile(
        {"562": 3, "9606": 1, "0": 2}, tmp_path / "gold.profile", "s1", "genus", total=6
    )
    text = path.read_text()
    assert "@SampleID:s1" in text
    assert "@Version:0.10.0" in text
    assert "562" in text and "9606" in text
    assert "\t0\t" not in text
    sections = parse_profile(text)
    assert sections[0][1][-1] == "PERCENTAGE"
    by_id = {r["TAXID"]: r for r in sections[0][2]}
    assert abs(float(by_id["562"]["PERCENTAGE"]) - 50.0) < 1e-6
