"""CAMISIM generate backend: configs, community design, hybrid read_type tags."""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

from samovar.camisim import (
    add_host_abundance,
    camisim_sizes_by_type,
    collect_generate_genomes,
    default_generate_config,
    design_communities,
    extract_read_type,
    filter_wgsim_zero_read_taxa,
    normalize_camisim_mode,
    run_from_config,
    setup_camisim_generate,
    tag_fastq_file,
    tag_read_id,
    types_for_mode,
    write_camisim_configs,
    write_distribution_files,
)
from samovar.combine_tables import combine_with_cpp, ensure_combine_binary
from samovar.parse_annotators import extract_read_type as parse_extract_read_type


REPO = Path(__file__).resolve().parents[1]
META = REPO / "data" / "test_genomes" / "meta"
HOST = REPO / "data" / "test_genomes" / "host" / "9606.fna"


def test_mode_aliases():
    assert normalize_camisim_mode("table") == "table"
    assert normalize_camisim_mode("illumina") == "illumina"
    assert normalize_camisim_mode("art") == "illumina"
    assert normalize_camisim_mode("ont") == "ont"
    assert normalize_camisim_mode("nanosim3") == "ont"
    assert normalize_camisim_mode("hybrid") == "hybrid"
    assert types_for_mode("illumina") == ["art"]
    assert types_for_mode("ont") == ["nanosim3"]
    assert types_for_mode("wgsim") == ["wgsim"]
    hybrid = types_for_mode("hybrid")
    assert "art" in hybrid and "nanosim3" in hybrid


def test_extract_and_tag_read_type():
    tagged = tag_read_id("@r1 taxid:562", "illumina", "562")
    assert "read_type:illumina" in tagged
    assert "taxid:562" in tagged
    assert extract_read_type(tagged) == "illumina"
    assert parse_extract_read_type("x|read_type:ont|taxid:9606") == "ont"
    assert extract_read_type("plain") == ""


def test_write_editable_camisim_config(tmp_path):
    args = type(
        "Args",
        (),
        {
            "genome_dir": str(META),
            "host_genome": str(HOST),
            "output_dir": str(tmp_path / "out"),
            "n_samples": 2,
            "total_reads": 400,
            "host_fraction": "0.1",
            "seed": 7,
            "model": "hiseq",
            "cores": 1,
            "camisim_mode": "table",
            "camisim_config": None,
            "size_gbp": None,
            "simulator": "camisim",
        },
    )()
    result = setup_camisim_generate(args)
    yaml_path = Path(result["config"])
    pipeline = Path(result["pipeline"])
    assert yaml_path.is_file()
    assert pipeline.is_file()
    cfg = yaml.safe_load(yaml_path.read_text())
    assert cfg["simulator"] == "camisim"
    assert cfg["mode"] == "table"
    assert cfg["n_samples"] == 2
    assert cfg["total_reads"] == 400
    assert cfg["seed"] == 7
    assert cfg["just_community_design"] is True
    assert Path(cfg["genome_dir"]).is_dir()
    text = pipeline.read_text()
    assert "samovar.camisim" in text
    assert "camisim.yaml" in text
    loc = Path(result["genome_locations"])
    assert loc.is_file()
    body = loc.read_text()
    assert "562" in body
    assert "9606" in body


def test_map_ncbi_id_for_old_taxdump(tmp_path):
    import tarfile

    from samovar.camisim import load_taxdump_taxids, map_ncbi_id_for_taxdump, write_genome_tables

    dump = tmp_path / "taxdump.tar.gz"
    nodes = "1\t|\t1\t|\tno rank\t|\n10847\t|\t1\t|\tspecies\t|\n562\t|\t1\t|\tspecies\t|\n"
    with tarfile.open(dump, "w:gz") as tar:
        info = tarfile.TarInfo("nodes.dmp")
        payload = nodes.encode()
        info.size = len(payload)
        tar.addfile(info, fileobj=__import__("io").BytesIO(payload))
    ids = load_taxdump_taxids(str(dump))
    assert "562" in ids
    assert map_ncbi_id_for_taxdump("562", ids) == "562"
    assert map_ncbi_id_for_taxdump("2886930", ids) == "10847"
    assert map_ncbi_id_for_taxdump("99999999", ids) == "1"
    rows = [
        {
            "genome_ID": "2886930",
            "path": "/tmp/x.fna",
            "NCBI_ID": "2886930",
            "OTU": "2886930",
            "novelty_category": "known_strain",
        }
    ]
    loc, meta = write_genome_tables(tmp_path / "cam", rows, taxdump=str(dump))
    assert "2886930\t2886930\t10847\tknown_strain" in meta.read_text()
    assert loc.is_file()


def test_sanitize_fasta_headers_for_wgsim(tmp_path):
    from samovar.camisim import sanitize_fasta_for_camisim

    src = tmp_path / "562.fna"
    src.write_text(">Ecoli.fna|taxid:562|-NC_000913.3 desc\nACGT\n")
    dest = tmp_path / "clean.fna"
    sanitize_fasta_for_camisim(src, dest, "562")
    text = dest.read_text()
    assert text.startswith(">562c1\n")
    assert "_" not in text.splitlines()[0]
    assert "ACGT" in text


def test_illumina_config_skips_iss_flag(tmp_path):
    cfg = default_generate_config(
        genome_dir=str(META),
        output_dir=str(tmp_path),
        host_genome=str(HOST),
        n_samples=2,
        total_reads=200,
        mode="illumina",
    )
    assert cfg["mode"] == "illumina"
    assert cfg["camisim_types"] == ["art"]
    assert cfg["just_community_design"] is False
    # 200 ART pairs * 2 ends * 150 bp, not the old 0.01-Gbp floor
    # that accidentally generated ~33,000 pairs.
    assert cfg["size_gbp"] == 0.00006
    assert cfg["size_gbp_by_type"] == {"art": 0.00006}
    paths = write_camisim_configs(cfg, str(tmp_path))
    nf = Path(paths["nextflow"]).read_text()
    assert 'type = "art"' in nf
    assert "just_community_design = false" in nf
    assert "distribution_files =" in nf
    assert "ncbi_taxdump_file" in nf


def test_community_design_reproducible(tmp_path):
    rows = collect_generate_genomes(str(META), str(HOST))
    ids = [r["genome_ID"] for r in rows]
    a = design_communities(ids, n_samples=3, seed=42)
    b = design_communities(ids, n_samples=3, seed=42)
    c = design_communities(ids, n_samples=3, seed=1)
    assert a == b
    assert a != c
    assert len(a) == 3
    for table in a:
        assert abs(sum(table.values()) - 1.0) < 1e-6
    files = write_distribution_files(tmp_path / "d", a)
    assert len(files) == 3
    first = files[0].read_text()
    assert "\t" in first
    assert not first.lower().startswith("genome")
    gid, abund = first.splitlines()[0].split("\t")
    assert gid in ids
    float(abund)


def test_type_sizes_split_hybrid_records_and_host_fraction():
    sizes = camisim_sizes_by_type(
        200,
        ["art", "nanosim3"],
        art_read_length=150,
        nanosim_read_length=4500,
    )
    assert sizes == {"art": 0.00003, "nanosim3": 0.00045}

    rows = [
        {"genome_ID": "562", "host": "0"},
        {"genome_ID": "9606", "host": "1"},
    ]
    samples = add_host_abundance([{"562": 1.0}], rows, "0.25")
    assert samples == [{"562": 0.75, "9606": 0.25}]


def test_wgsim_filters_taxa_that_would_receive_zero_pairs(tmp_path):
    tiny = tmp_path / "tiny.fna"
    large = tmp_path / "large.fna"
    tiny.write_text(">tiny\n" + "A" * 10 + "\n")
    large.write_text(">large\n" + "A" * 10000 + "\n")
    rows = [
        {"genome_ID": "tiny", "path": str(tiny)},
        {"genome_ID": "large", "path": str(large)},
    ]
    filtered = filter_wgsim_zero_read_taxa(
        [{"tiny": 0.01, "large": 0.99}],
        rows,
        size_gbp=0.00003,
        read_length=150,
    )
    assert filtered == [{"large": 1.0}]


def test_hybrid_read_type_survives_combiner(tmp_path):
    ensure_combine_binary()
    reports = tmp_path / "reports"
    reports.mkdir()
    ill = "readA|taxid:562|read_type:illumina"
    ont = "readB|taxid:9606|read_type:ont"
    (reports / "1_kaiju.kaiju.out").write_text(f"C\t{ill}\t562\nC\t{ont}\t9606\n")
    (reports / "1_kraken2.kraken2.out").write_text(
        f"C\t{ill}\torg (taxid 562)\t100|100\t0:1\n"
        f"C\t{ont}\torg (taxid 9606)\t80|80\t0:1\n"
    )
    out = tmp_path / "ann"
    combine_with_cpp(str(reports), str(out), split_n=1)
    csv = next(out.glob("*.annotation.csv"))
    df = pd.read_csv(csv)
    assert "read_type" in df.columns
    by = df.set_index("seq")
    assert str(by.loc[ill, "read_type"]).lower() == "illumina"
    assert str(by.loc[ont, "read_type"]).lower() == "ont"


def test_tag_fastq_file(tmp_path):
    src = tmp_path / "in.fastq"
    src.write_text("@r1 extra\nACGT\n+\nIIII\n")
    dest = tmp_path / "out.fastq"
    tag_fastq_file(src, dest, "ont", "4932")
    text = dest.read_text()
    assert "read_type:ont" in text
    assert "taxid:4932" in text
    assert text.splitlines()[1] == "ACGT"
    tag_fastq_file(dest, dest, "ont", "4932")
    assert dest.read_text().splitlines()[1] == "ACGT"


def test_score_annotators_splits_hybrid_read_type():
    from samovar.scores import score_annotators

    work = pd.DataFrame(
        {
            "kaiju": ["562"] * 6 + ["9606"] * 6,
            "true": ["562"] * 6 + ["9606"] * 6,
            "read_type": ["illumina"] * 3 + ["ont"] * 3 + ["illumina"] * 3 + ["ont"] * 3,
        }
    )
    table = score_annotators(work, ["kaiju"])
    assert "read_type" in table.columns
    types = set(table["read_type"].astype(str))
    assert "all" in types
    assert "illumina" in types
    assert "ont" in types


def test_table_mode_iss_from_distributions(tmp_path, monkeypatch):
    from samovar import camisim as cam

    calls = []

    def fake_genome(*args, **kwargs):
        prefix = Path(args[1] if args else kwargs.get("output_file"))
        prefix.parent.mkdir(parents=True, exist_ok=True)
        (Path(str(prefix) + "_R1.fastq")).write_text("@h\nA\n+\nI\n")
        (Path(str(prefix) + "_R2.fastq")).write_text("@h\nA\n+\nI\n")
        calls.append("host")
        return str(prefix) + "_R1.fastq", str(prefix) + "_R2.fastq"

    def fake_meta(**kwargs):
        out = Path(kwargs["output_dir"])
        (out / "meta_full_R1.fastq").write_text("@m|taxid:562\nACGT\n+\nIIII\n")
        (out / "meta_full_R2.fastq").write_text("@m|taxid:562\nACGT\n+\nIIII\n")
        calls.append("meta")
        return str(out / "meta_full_R1.fastq"), str(out / "meta_full_R2.fastq")

    monkeypatch.setattr("samovar.table2iss.generate_reads_genome", fake_genome)
    monkeypatch.setattr("samovar.table2iss.generate_reads_metagenome", fake_meta)

    args = type(
        "Args",
        (),
        {
            "genome_dir": str(META),
            "host_genome": str(HOST),
            "output_dir": str(tmp_path / "run"),
            "n_samples": 2,
            "total_reads": 20,
            "host_fraction": "0.2",
            "seed": 1,
            "model": "hiseq",
            "cores": 1,
            "camisim_mode": "table",
            "camisim_config": None,
            "size_gbp": 0.0001,
            "simulator": "camisim",
        },
    )()
    setup = setup_camisim_generate(args)
    result = run_from_config(setup["config"])
    assert result["mode"] == "table"
    assert result["used_iss"] is True
    assert Path(result["abundance"]).is_file()
    reads = [Path(p) for p in result["reads"]]
    assert any(p.name.endswith("_full_R1.fastq") for p in reads)
    assert (tmp_path / "run" / "initial" / "1_full_R1.fastq").is_file()
    assert calls


def test_camisim_config_overlay_keeps_nested_template(tmp_path):
    template = tmp_path / "user.yaml"
    template.write_text(
        "distribution:\n  log_mu: 3.5\n  log_sigma: 0.5\nart:\n  fragment_size_mean: 400\ncustom_flag: extra\n",
        encoding="utf-8",
    )
    args = type(
        "Args",
        (),
        {
            "genome_dir": str(META),
            "host_genome": str(HOST),
            "output_dir": str(tmp_path / "out"),
            "n_samples": 2,
            "total_reads": 400,
            "host_fraction": "0.1",
            "seed": 7,
            "model": "hiseq",
            "cores": 1,
            "camisim_mode": "illumina",
            "camisim_config": str(template),
            "size_gbp": None,
            "simulator": "camisim",
        },
    )()
    result = setup_camisim_generate(args)
    cfg = yaml.safe_load(Path(result["config"]).read_text())
    assert cfg["n_samples"] == 2
    assert cfg["distribution"]["log_mu"] == 3.5
    assert cfg["distribution"]["log_sigma"] == 0.5
    assert cfg["distribution"]["mode"] == "differential"
    assert cfg["art"]["fragment_size_mean"] == 400
    assert cfg["custom_flag"] == "extra"


def test_table_nextflow_config_leaves_distributions_empty(tmp_path):
    cfg = default_generate_config(
        genome_dir=str(META),
        output_dir=str(tmp_path),
        host_genome=str(HOST),
        n_samples=2,
        total_reads=200,
        mode="table",
    )
    paths = write_camisim_configs(cfg, str(tmp_path))
    nf = Path(paths["nextflow"]).read_text()
    assert "just_community_design = true" in nf
    assert 'distribution_files = ""' in nf


def test_nextflow_run_cmd_uses_exclusive_config(tmp_path):
    from samovar.camisim import nextflow_run_cmd

    nf = tmp_path / "nextflow.config"
    nf.write_text("//\n")
    cmd = nextflow_run_cmd("nextflow", "/opt/CAMISIM", nf, tmp_path / "work")
    assert cmd[:3] == ["nextflow", "run", "main.nf"]
    assert "--config" in cmd
    assert "-c" not in cmd
    assert str(nf.resolve()) in cmd


def test_nextflow_config_uses_nanosim_sidecar(tmp_path, monkeypatch):
    from samovar.camisim import default_generate_config, write_camisim_configs

    env = tmp_path / "ns"
    (env / "bin").mkdir(parents=True)
    exe = env / "bin" / "simulator.py"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    monkeypatch.setattr(
        "samovar.camisim.tool_env_prefix",
        lambda name, cfg=None: str(env) if "nano" in name else None,
    )
    monkeypatch.setattr("samovar.camisim.discover_nanosim", lambda: str(exe))
    cfg = default_generate_config(
        genome_dir=str(META),
        output_dir=str(tmp_path),
        host_genome=str(HOST),
        n_samples=1,
        total_reads=200,
        mode="ont",
    )
    paths = write_camisim_configs(cfg, str(tmp_path))
    nf = Path(paths["nextflow"]).read_text()
    assert 'type = "nanosim3"' in nf
    assert env.as_posix() in nf.replace("\\", "/")
    assert "simulate_reads_.*nanosim3" in nf


def test_harvest_prefers_pooled_mates(tmp_path):
    from samovar.camisim import _harvest_camisim_reads, _mate_kind

    raw = tmp_path / "raw" / "camisim_raw" / "sample_0" / "reads" / "fastq"
    raw.mkdir(parents=True)
    (raw / "sample_0_01.fq").write_text("@p1\nAC\n+\nII\n")
    (raw / "sample_0_02.fq").write_text("@p2\nGT\n+\nII\n")
    (raw / "sample0_5621.fq").write_text("@g1\nAA\n+\nII\n")
    (raw / "sample0_5622.fq").write_text("@g2\nTT\n+\nII\n")
    rows = [{"genome_ID": "562", "NCBI_ID": "562"}]
    assert _mate_kind(raw / "sample_0_02.fq", rows) == "R2"
    out = tmp_path / "initial"
    written = _harvest_camisim_reads(tmp_path / "raw", out, rows, "wgsim", 1)
    r1 = (out / "1_wgsim_R1.fastq").read_text()
    r2 = (out / "1_wgsim_R2.fastq").read_text()
    assert "p1" in r1 and "g1" not in r1
    assert "p2" in r2 and "g2" not in r2
    assert written


def test_promote_removes_tech_fastq_copies(tmp_path):
    from samovar.camisim import _promote_to_full

    out = tmp_path / "initial"
    out.mkdir()
    src = out / "1_wgsim_R1.fastq"
    src.write_text("@p1|read_type:wgsim\nAC\n+\nII\n")
    (out / "1_wgsim_R2.fastq").write_text("@p2|read_type:wgsim\nGT\n+\nII\n")
    written = _promote_to_full([src], out, 1, "wgsim")
    assert (out / "1_full_R1.fastq").is_file()
    assert "read_type:wgsim" in (out / "1_full_R1.fastq").read_text()
    assert not (out / "1_wgsim_R1.fastq").exists()
    assert not (out / "1_wgsim_R2.fastq").exists()
    assert any(Path(p).name == "1_full_R1.fastq" for p in written)


def test_hybrid_flattens_paired_and_single_end_reads(tmp_path):
    from samovar.camisim import _merge_hybrid_reads

    out = tmp_path / "initial"
    out.mkdir()
    (out / "1_art_R1.fastq").write_text("@a1\nAC\n+\nII\n")
    (out / "1_art_R2.fastq").write_text("@a2\nGT\n+\nII\n")
    (out / "1_nanosim3_R1.fastq").write_text("@n1\nAAAA\n+\nIIII\n")
    (out / "1_nanosim3_R2.fastq").write_text("")

    _merge_hybrid_reads({"art": [], "nanosim3": []}, out, 1)

    full_r1 = (out / "1_full_R1.fastq").read_text()
    assert "@a1" in full_r1 and "@a2" in full_r1 and "@n1" in full_r1
    assert (out / "1_full_R2.fastq").read_text() == ""

