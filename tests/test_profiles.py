"""Independent sample profiles, GT semantics and correction through public contracts."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from samovar.abundance import input_to_abundance_tables, load_table_input
from samovar.annotation_convert import dump_builtin, load_builtin
from samovar.parse_annotators import Annotation
from samovar.profiles import (change_basis, load_profile, metadata, read_counts,
                              read_profiles, validate, write_profiles)
from samovar.scores import profile_scores


def table(values, ids=("562", "4932"), sample="s"):
    return pd.DataFrame({"taxid": list(ids), "N_" + sample: list(values)})


@pytest.fixture
def external_kraken_profiler(tmp_path, monkeypatch):
    import os
    import shutil
    import sys

    folder = tmp_path / "external profiler space's Ω"
    folder.mkdir()
    wrapper = folder / "kraken2_profile.py"
    shutil.copyfile(Path(__file__).parent / "data/kraken2_profile_wrapper.py", wrapper)
    fake = folder / "native kraken2"
    fake.write_text(f'#!{sys.executable}\n' + '''import json, sys
from pathlib import Path
args = sys.argv[1:]
db = Path(args[args.index('--db') + 1])
(db / 'argv.json').write_text(json.dumps(args))
if (db / 'fail').exists():
    raise SystemExit(7)
Path(args[args.index('--report') + 1]).write_text((db / 'report').read_text())
''')
    fake.chmod(0o755)
    config = Path(os.environ["SAMOVAR_CONFIG"])
    cfg = json.loads(config.read_text())
    cfg['tools'] = {key: value for key, value in cfg['tools'].items() if key.split(':')[0] != 'kraken2'}
    cfg.get('tool_envs', {}).pop('kraken2', None)
    cfg["tools"]["kraken2"] = {"type": "annotator", "exec": {"path": str(fake)}}
    config.write_text(json.dumps(cfg))
    foreign = folder / "foreign"
    foreign.mkdir()
    for name in ("python", "python3", "kraken2"):
        path = foreign / name
        path.write_text("#!/bin/sh\nexit 91\n")
        path.chmod(0o755)
    monkeypatch.setenv("PATH", str(foreign) + os.pathsep + os.environ["PATH"])
    db = folder / "database"
    db.mkdir()
    (db / "report").write_text("30\t3\t3\tU\t0\tunclassified\n"
                              "70\t7\t0\tR\t1\troot\n"
                              "70\t7\t3\tG\t561\t  genus\n"
                              "40\t4\t4\tS\t562\t    species\n")
    reads = folder / "reads.fastq"
    reads.write_text("@read\nACGT\n+\nIIII\n")
    return wrapper, db, reads


@pytest.mark.parametrize("paired", [False, True])
def test_external_profiler_import_and_production_launch(external_kraken_profiler, paired, capsys):
    import subprocess
    from samovar.tools_import import main as import_main
    from samovar.annotators_wrapper import get_annotator_instance
    from argparse import Namespace
    from samovar.config import PipelineConfig

    wrapper, db, reads = external_kraken_profiler
    argv = ["--name", "external_profile", "--type", "profiler", "--exec-path", str(wrapper),
            "--tool-version", "test", "--db-path", str(db)]
    if paired:
        argv.append("--pytest")
    assert import_main(argv) == 0
    settings = wrapper.parent / "prepare.yaml"
    settings.write_text(__import__("yaml").safe_dump({"annotators": [dict(
        type="external_profile", cmd="external_profile", run_name="named_profile", db_path=str(db))]}))
    cfg = PipelineConfig.from_args(Namespace(input_config=str(settings), input_dir=None,
                                            output_dir=str(wrapper.parent / "run"), kraken2=None, kaiju=None))
    ann = cfg.annotators[0]
    assert ann.output_level == "profile"
    assert "not on PATH" not in capsys.readouterr().err
    tool = get_annotator_instance(ann.type, vars(ann), {})
    out = Path(tool.get_expected_outputs("patient_A", str(wrapper.parent))[0])
    command = tool.get_snakemake_shell_cmd(str(reads), str(reads) if paired else "", [str(out)])
    subprocess.run(["bash", "-e", "-c", command], check=True)
    frame, info = load_profile(out)
    assert dict(zip(frame.taxid, frame.N_patient_A)) == {"0": 3, "561": 3, "562": 4}
    assert info["identity"] == "named_profile" and info["database"] == str(db)
    assert info["rank"] == "taxid" and info["normalization"] == "all"
    args = json.loads((db / "argv.json").read_text())
    assert ("--paired" in args) == paired
    assert not list((wrapper.parent / ".tmp").iterdir())


@pytest.mark.parametrize("mode", ["empty", "no_hits", "blank", "malformed", "negative", "duplicate", "failure"])
def test_external_profiler_empty_and_failed_results(external_kraken_profiler, mode):
    import subprocess
    from samovar.paths import python_path
    wrapper, db, reads = external_kraken_profiler
    if mode == "empty":
        reads.write_text("")
    elif mode == "no_hits":
        (db / "report").write_text("100\t1\t1\tU\t0\tunclassified\n")
    elif mode == "failure":
        (db / "fail").touch()
    else:
        report = {"blank": "", "malformed": "broken\n", "negative": "0\t1\t-1\tS\t562\tspecies\n",
                  "duplicate": "50\t1\t1\tS\t562\tspecies\n" * 2}[mode]
        (db / "report").write_text(report)
    out = wrapper.parent / "result.csv"
    result = subprocess.run([python_path(), str(wrapper), "-i", str(reads), "-d", str(db),
                             "-o", str(out), "-t", "2"], capture_output=True, text=True)
    if mode in {"empty", "no_hits"}:
        assert result.returncode == 0, result.stderr
        frame, _ = load_profile(out)
        assert frame.empty if mode == "empty" else frame.taxid.tolist() == ["0"]
    else:
        assert result.returncode != 0
        assert not out.exists() and not Path(str(out) + ".profile.json").exists()
    assert not list((wrapper.parent / ".tmp").iterdir())


def test_sylph_native_profile_uses_explicit_taxonomy_and_organism_abundance(tmp_path):
    mapping = tmp_path / "genomes.syldb.taxonomy.tsv"
    mapping.write_text("Genome_file\ttaxid\ng1.fna\t562\ng2.fna\t562\n")
    raw = tmp_path / "native.tsv"
    header = "Sample_file\tGenome_file\tTaxonomic_abundance\tSequence_abundance\n"
    raw.write_text(header + "reads_R1.fq\tg1.fna\t25\t60\nreads_R1.fq\tg2.fna\t75\t40\n")
    settings = {"genome_taxonomy": str(mapping), "sample": "original_sample"}
    frame, info = load_profile(raw, settings)
    assert frame.to_dict("list") == {"taxid": ["562"], "N_original_sample": [100.]}
    assert (info["unit"], info["basis"], info["normalization"]) == ("percent", "organisms", "classified")
    assert info["rank"] == "taxid"
    mapping.write_text("Genome_file\ttaxid\ng1.fna\t511145\ng2.fna\t83333\n")
    from samovar.profiles import at_rank
    frame, info = load_profile(raw, settings)
    species, _ = at_rank(frame, info, "species", Path(__file__).parent / "data/taxdump")
    assert species.taxid.tolist() == ["562"] and species.N_original_sample.item() == 100
    mapping.write_text("Genome_file\ttaxid\ng1.fna\t562\ng2.fna\t562\n")
    raw.write_text(header)
    frame, _ = load_profile(raw, {**settings, "format": "sylph"})
    assert frame.empty and list(frame) == ["taxid", "N_original_sample"]
    raw.write_text(header + "reads_R1.fq\tunknown.fna\t100\t100\n")
    with pytest.raises(ValueError, match="lack explicit taxonomy"):
        load_profile(raw, settings)
    raw.write_text(header + "a\tg1.fna\t50\t50\nb\tg2.fna\t50\t50\n")
    with pytest.raises(ValueError, match="one Sylph sample"):
        load_profile(raw, settings)
    raw.write_text(header + "a\tg1.fna\t-1\t100\n")
    with pytest.raises(ValueError, match="nonnegative"):
        load_profile(raw, settings)
    raw.write_text(header + "a\tg1.fna\t33.3334\t0\na\tg2.fna\t33.3334\t0\na\tg1.fna\t33.3333\t0\n")
    frame, _ = load_profile(raw, settings)
    assert frame.N_original_sample.sum() == pytest.approx(100)
    with pytest.raises(ValueError, match="NCBI abundance"):
        load_profile(raw, {**settings, "taxonomy": "gtdb"})
    mapping.write_text("Genome_file\ttaxid\ng1.fna\tSGB1\ng2.fna\t562\n")
    with pytest.raises(ValueError, match="positive NCBI taxids"):
        load_profile(raw, settings)


def annotation(frame, name="mpa", **info):
    return Annotation.from_abundance_tables({name: frame}, {name: metadata(info)})


def test_mixed_annotation_keeps_both_channels_and_roundtrip(tmp_path):
    obj = Annotation.from_long_table(pd.DataFrame({"seq": ["r1", "r2"], "sample": ["s", "s"], "taxID_kaiju_0": [562, 562]}))
    obj.add_profiles({"mpa_ref": table([25, 75])}, {"mpa_ref": metadata(dict(unit="percent", basis="organisms", database="ref"))})
    got = input_to_abundance_tables(obj)
    assert set(got) == {"kaiju", "mpa_ref"}
    assert len(obj.DataFrame) == 2 and "seq" not in got["mpa_ref"]
    dest = tmp_path / "profiles"
    dump_builtin(obj, dest, "abundance")
    loaded = load_table_input(dest)
    assert loaded.DataFrame.empty
    assert loaded.profile_metadata["mpa_ref"]["basis"] == "organisms"
    pd.testing.assert_frame_equal(input_to_abundance_tables(loaded)["mpa_ref"], got["mpa_ref"], check_dtype=False)
    # Conversion caches must not introduce a second copy of the read output.
    assert set(input_to_abundance_tables(obj)) == set(got)


@pytest.mark.parametrize("stage", ["initial", "regenerated", "reprofiled"])
def test_report_contract_preserves_mixed_profiles_without_truth(tmp_path, stage):
    from samovar.profiles import report_profiles
    obj = Annotation.from_long_table(pd.DataFrame({"seq": ["r1", "r2"], "sample": ["s", "s"],
        "taxID_kaiju_0": [562, 0], "taxID_kaiju_1": [4932, 4932]}))
    obj.add_profiles({"mpa": table([20, 70])}, {"mpa": metadata(dict(unit="percent", basis="organisms",
        database="markers", unassigned_by_sample={"s": 0.1}))})
    obj.abundance_tables["kaiju"] = table([2, 0])  # old conversion cache, not a new profiler
    config = {"run_config": [dict(run_name="kaiju", type="kaiju", db_path="refseq")]}
    rows = report_profiles(obj, tmp_path, stage, config)
    tables, metas = read_profiles(tmp_path)
    assert set(tables) == {"taxID_kaiju_0", "taxID_kaiju_1", "mpa"}
    assert tables["taxID_kaiju_0"].set_index("taxid").N_s.to_dict() == {"0": 1, "562": 1}
    assert metas["taxID_kaiju_0"]["unassigned_by_sample"] == {"s": 0.5}
    assert metas["taxID_kaiju_1"]["database"] == "refseq"
    assert metas["mpa"]["unit"] == "percent" and metas["mpa"]["basis"] == "organisms"
    assert all(row["stage"] == stage for row in rows.values())
    assert len(rows) == 3 and "seq" not in tables["mpa"]
    # Independent units/bases stay separate, and rebuilding cannot change counts.
    before = {p.name: p.read_bytes() for p in tmp_path.iterdir()}
    report_profiles(obj, tmp_path, stage, config)
    assert before == {p.name: p.read_bytes() for p in tmp_path.iterdir()}


def test_report_contract_empty_unclassified_unresolved_and_invalid_metadata(tmp_path):
    from samovar.profiles import report_profiles
    obj = Annotation.from_abundance_tables(
        {"empty": table([], ids=()), "unknown": table([5], ids=("0",)),
         "native": table([100], ids=("unmapped:SGB1",))},
        {"empty": metadata(), "unknown": metadata({"normalization": "all"}),
         "native": metadata({"unit": "percent", "basis": "organisms"})})
    rows = list(report_profiles(obj, tmp_path, "initial").values())
    assert [r["status"] for r in rows] == ["empty", "all-unclassified", "available"]
    assert rows[2]["unresolved_fraction"] == 1
    assert rows[0]["unassigned_fraction"] is None
    obj.profile_metadata["empty"]["unassigned_by_sample"] = {"s": 1.0}
    assert report_profiles(obj, tmp_path, "initial")["0"]["status"] == "all-unclassified"
    obj.profile_metadata["native"]["unassigned_by_sample"] = {"absent": 0.1}
    with pytest.raises(ValueError, match="unknown samples"):
        report_profiles(obj, tmp_path, "initial")


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -1, "broken"])
def test_invalid_values_are_errors(value):
    with pytest.raises((ValueError, TypeError)):
        validate(table([value, 1]))


def test_presence_is_binary_and_empty_sample_is_preserved(tmp_path):
    with pytest.raises(ValueError, match="0 or 1"):
        validate(table([0.1, 1]), {"kind": "presence"})
    frame = table([], ids=())
    write_profiles(tmp_path, {"gt": frame}, {"gt": metadata({"kind": "presence"})})
    tables, metas = read_profiles(tmp_path)
    assert list(tables["gt"]) == ["taxid", "N_s"] and tables["gt"].empty
    scores = profile_scores(table([1, 0]), tables["gt"], metadata(), metas["gt"]).iloc[0]
    assert scores.fp == 1 and scores.fn == scores.tp == 0
    assert np.isnan(scores.recall) and np.isnan(scores.l1_norm)


def test_profile_scores_count_taxa_not_artificial_reads():
    pred = table([1, 99], ids=("562", "9606"))
    gt = table([1, 1], ids=("562", "4932"))
    score = profile_scores(pred, gt, metadata(), metadata({"kind": "presence"})).iloc[0]
    assert (score.tp, score.fp, score.fn) == (1, 1, 1)
    assert score.precision == score.recall == score.f1 == 0.5
    assert np.isnan(score.bray_curtis)
    with pytest.raises(ValueError, match="sample identities"):
        profile_scores(pred, gt.rename(columns={"N_s": "N_other"}), metadata(), metadata())
    with pytest.raises(ValueError, match="basis"):
        profile_scores(pred, gt, metadata(), metadata({"basis": "organisms"}))


def test_length_conversion_and_exact_fragment_allocation():
    frame = table([50, 50])
    info = metadata(dict(unit="percent", basis="organisms"))
    reads, rm = change_basis(frame, info, "reads", {"562": 1000, "4932": 3000})
    assert reads.N_s.tolist() == [0.25, 0.75]
    back, _ = change_basis(reads, rm, "organisms", {"562": 1000, "4932": 3000})
    assert back.N_s.tolist() == [0.5, 0.5]
    counts, cm = read_counts(frame, info, 101, {"562": 1000, "4932": 3000})
    assert counts.N_s.tolist() == [25, 76] and cm["basis"] == "reads"
    with pytest.raises(ValueError, match="reference lengths"):
        read_counts(frame, info, 100, {"562": 1000})


def test_profile_reference_selection_is_pinned_and_changes_fail(tmp_path, monkeypatch):
    from samovar.profiles import reference_lengths
    from samovar.table2iss import _resolve_genomes_for_taxids
    source = tmp_path / "reference.fa"
    source.write_text(">chromosome\nACGTACGT\n")
    monkeypatch.setattr("samovar.genome_resolve.resolve_genome_file", lambda *a, **kw: str(source))
    selected = tmp_path / "genomes"
    assert reference_lengths({"562"}, selected) == {"562": 8}
    monkeypatch.setattr("samovar.genome_resolve.resolve_genome_file", lambda *a, **kw: pytest.fail("must reuse selected reference"))
    assert _resolve_genomes_for_taxids(["562"], str(selected), "", True) == {"562": str(source)}
    source.write_text(">chromosome\nACGT\n")
    with pytest.raises(ValueError, match="reference changed"):
        _resolve_genomes_for_taxids(["562"], str(selected), "", True)


def test_incompatible_taxonomy_and_unassigned_length_are_explicit_errors():
    with pytest.raises(ValueError, match="taxonomy"):
        profile_scores(table([1, 2]), table([1, 2]), metadata(), metadata({"taxonomy": "gtdb"}))
    with pytest.raises(ValueError, match="unassigned mass"):
        change_basis(table([30, 50]), metadata(dict(unit="percent", basis="organisms", normalization="all")), "reads", {"562": 1000, "4932": 3000})


def test_native_metaphlan_retains_unresolved_taxa_and_avoids_clade_doublecount(tmp_path):
    path = tmp_path / "mpa.txt"
    path.write_text("#clade_name\tNCBI_tax_id\trelative_abundance\n"
                    "k__Bacteria\t2\t100\n"
                    "k__Bacteria|s__E_coli\t2|562\t70\n"
                    "k__Bacteria|s__unknown\t2|\t30\n")
    frame, meta = load_profile(path, {"sample": "real_sample"})
    assert meta["basis"] == "organisms" and meta["unit"] == "percent"
    assert frame.N_real_sample.sum() == 100 and len(frame) == 2
    assert any(t.startswith("unmapped:") for t in frame.taxid)


def test_cami_multiple_samples_including_empty(tmp_path):
    path = tmp_path / "gold.profile"
    path.write_text("@SampleID:empty\n@TaxonomyID:ncbi\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n"
                    "@SampleID:other\n@@TAXID\tRANK\tTAXPATH\tTAXPATHSN\tPERCENTAGE\n"
                    "562\tspecies\t562\tE_coli\t100\n")
    frame, meta = load_profile(path)
    assert list(frame) == ["taxid", "N_empty", "N_other"]
    assert frame.N_empty.sum() == 0 and frame.N_other.sum() == 100
    assert load_builtin(path, "cami").DataFrame.empty


def test_typed_cami_roundtrip_preserves_identity_units_and_unassigned_mass(tmp_path):
    obj = annotation(table([30, 50]), unit="percent", basis="organisms", normalization="all", database="marker-db")
    path = dump_builtin(obj, tmp_path / "export.profile", "cami")
    back = load_builtin(path, "cami")
    assert set(back.abundance_tables) == {"mpa"}
    assert back.profile_metadata["mpa"]["basis"] == "organisms"
    assert back.profile_metadata["mpa"]["database"] == "marker-db"
    assert back.abundance_tables["mpa"].N_s.sum() == 80


@pytest.mark.parametrize("mode, expected", [("table", [0.25, 0.75]), ("illumina", [0.5, 0.5])])
def test_camisim_boundary_does_not_apply_genome_length_twice(tmp_path, monkeypatch, mode, expected):
    from samovar import camisim
    from samovar.reads_generators import _simulate_camisim_sample
    import yaml
    genomes = {}
    for tax, length in (("562", 1000), ("4932", 3000)):
        path = tmp_path / (tax + ".fna")
        path.write_text(">chr\n" + "A" * length + "\n")
        genomes[tax] = str(path)
    monkeypatch.setattr(camisim, "default_generate_config", lambda **kw: {"mode": mode})
    def run(path):
        cfg = yaml.safe_load(Path(path).read_text())
        assert [cfg["communities"][0][tax] for tax in ("562", "4932")] == expected
        # Stop at the simulator boundary after checking the actual serialized input.
        raise RuntimeError("boundary checked")
    monkeypatch.setattr(camisim, "run_from_config", run)
    with pytest.raises(RuntimeError, match="boundary checked"):
        _simulate_camisim_sample("camisim", {"562": 25, "4932": 75}, genomes, tmp_path / "out", "s", "mpa",
            model="hiseq", cpus=1, seed=42, gzip_reads=False, flags=[], ids=[], config={})


def test_ratio_correction_saved_apply_and_held_out_bias(tmp_path):
    from samovar.abundance_correctors import run_export
    train = annotation(table([80, 20, 10], ids=("562", "4932", "9606")))
    truth = table([40, 60])
    held_out = annotation(table([40, 10, 5], ids=("562", "4932", "9606")))
    target = table([20, 30])
    model = tmp_path / "model.json"
    cfg = dict(reference=train, reference_truth=(truth, metadata()), profile_model=str(model))
    run_export("profile_ratio", held_out, tmp_path / "first", cfg)
    corrected, metas = read_profiles(tmp_path / "first")
    result = corrected["mpa_profile_ratio"]
    before = profile_scores(held_out.abundance_tables["mpa"], target, metadata(), metadata()).iloc[0]
    after = profile_scores(result, target, metas["mpa_profile_ratio"], metadata()).iloc[0]
    assert after.l1_norm < before.l1_norm
    assert result.loc[result.taxid == "9606", "N_s"].item() == 0
    payload = json.loads(model.read_text())
    assert payload["profiles"]["mpa"]["factors"]["9606"] == 0
    run_export("profile_ratio", held_out, tmp_path / "loaded", {"profile_model": str(model)})
    pd.testing.assert_frame_equal(read_profiles(tmp_path / "loaded")[0]["mpa_profile_ratio"], result)
    unseen = annotation(table([40, 10, 5, 2], ids=("562", "4932", "9606", "2886930")))
    run_export("profile_ratio", unseen, tmp_path / "unseen", {"profile_model": str(model)})
    assert read_profiles(tmp_path / "unseen")[0]["mpa_profile_ratio"].query("taxid == '2886930'").N_s.item() > 0
    with pytest.raises(ValueError, match="basis"):
        run_export("profile_ratio", annotation(table([50, 50]), basis="organisms"), tmp_path / "bad", {"profile_model": str(model)})


def test_profile_input_prepare_does_not_create_per_read_model(tmp_path):
    from samovar.config import PipelineConfig, AnnotatorConfig
    c = PipelineConfig(output_dir=str(tmp_path), annotators=[AnnotatorConfig(
        run_name="mpa_ref", type="custom", cmd="unused", db_path=".", output_level="profile")])
    script = Path(c.generate_pipeline(str(tmp_path))).read_text()
    assert "samovar.profiles collect" in script and "samovar.profiles correct" in script
    assert "ML.py" not in script and "samovar.combine_tables finalize" not in script
    cfg = __import__("yaml").safe_load((tmp_path / ".log/configs/config_reannotate.yaml").read_text())
    assert cfg["run_config"][0]["output_level"] == "profile"
    imported = PipelineConfig(output_dir=str(tmp_path / "imported"), annotators=[],
                              profile_inputs={"external": {"path": "profile.csv"}})
    with pytest.raises(ValueError, match="matching run_name"):
        imported.generate_pipeline(imported.output_dir)
    imported.endpoint = "viz_initial"
    script = Path(imported.generate_pipeline(imported.output_dir)).read_text()
    assert "Using supplied sample profiles" in script

def test_bracken_profile_uses_estimated_read_counts(tmp_path):
    raw = tmp_path / 'sample.bracken.tsv'
    header = 'name\ttaxonomy_id\ttaxonomy_lvl\tnew_est_reads\tfraction_total_reads\n'
    raw.write_text(header + 'Ecoli\t562\tS\t25\t0.9\nYeast\t4932\tS\t75\t0.1\n')
    frame, info = load_profile(raw, {'sample': 'original'})
    assert frame.to_dict('list') == {'taxid': ['562', '4932'], 'N_original': [25., 75.]}
    assert (info['unit'], info['basis'], info['rank']) == ('counts', 'reads', 'species')
    raw.write_text(header)
    frame, _ = load_profile(raw, {'format': 'bracken', 'sample': 'original'})
    assert frame.empty and list(frame) == ['taxid', 'N_original']
    raw.write_text(header + 'Ecoli\t562\tS\t-1\t0.9\n')
    with pytest.raises(ValueError, match='nonnegative'):
        load_profile(raw, {'format': 'bracken'})
    raw.write_text(header + 'Ecoli\t562\tG\t1\t0.9\n')
    with pytest.raises(ValueError, match='species'):
        load_profile(raw, {'format': 'bracken'})


def test_bracken_estimator_configured_python_empty_threshold_and_failure(tmp_path, monkeypatch):
    import subprocess
    from samovar.profiles import run_bracken
    from samovar.paths import python_path
    report = tmp_path / 'sample report.txt'
    report.write_text('100\t12\t0\tR\t1\troot\n100\t12\t12\tS\t562\t  Ecoli\n')
    (tmp_path / 'database126mers.kmer_distrib').write_text('distribution')
    output = tmp_path / 'raw.tsv'
    script = str(tmp_path / 'program with spaces.py')
    Path(script).touch()
    monkeypatch.setattr('samovar.paths.resolve_executable', lambda *a, **k: script)
    calls = []
    def run(argv, **kwargs):
        calls.append(argv)
        assert kwargs['check']
        raise subprocess.CalledProcessError(2, argv)
    monkeypatch.setattr(subprocess, 'run', run)
    with pytest.raises(subprocess.CalledProcessError):
        run_bracken(report, tmp_path, output, 126, 10)
    assert calls[0][:2] == [python_path(), script]
    assert calls[0][calls[0].index('-i')+1] == str(report)
    calls.clear()
    run_bracken(report, tmp_path, output, 126, 13)
    assert not calls and load_profile(output, {'format': 'bracken'})[0].empty
    report.write_text('broken report\n')
    with pytest.raises(ValueError, match='Malformed'):
        run_bracken(report, tmp_path, output, 126)
    with pytest.raises(FileNotFoundError):
        run_bracken(report, tmp_path, output, 100)


def test_motus3_profile_preserves_species_unknowns_mass_and_roundtrip(tmp_path):
    raw = tmp_path / 'native.motus'
    raw.write_text('# motus version 3.1.0\n#consensus_taxonomy\tNCBI_tax_id\toriginal\n'
                   'Ecoli [ref_mOTU_v3_1]\t562\t0.2\n'
                   'Other Ecoli [ref_mOTU_v3_2]\t562\t0.1\n'
                   'unknown Clostridiales [meta_mOTU_v3_3]\t186802\t0.3\n'
                   'Novel [ext_mOTU_v3_4]\tNA\t0.1\nunassigned\tNA\t0.3\n')
    frame, info = load_profile(raw)
    assert frame.set_index('taxid').N_original.to_dict() == {
        '562': pytest.approx(0.3), 'unmapped:meta_mOTU_v3_3': 0.3,
        'unmapped:ext_mOTU_v3_4': 0.1, '0': 0.3}
    assert (info['unit'], info['basis'], info['normalization']) == ('fraction', 'organisms', 'all')
    assert info['unassigned_by_sample'] == {'original': 0.3}
    obj = Annotation.from_abundance_tables({'motus_db': frame}, {'motus_db': info})
    dump_builtin(obj, tmp_path/'exported', 'abundance')
    loaded = load_table_input(tmp_path/'exported')
    pd.testing.assert_frame_equal(input_to_abundance_tables(loaded)['motus_db'], frame, check_dtype=False)
    assert loaded.profile_metadata['motus_db']['native_header'] == info['native_header']
    raw.write_text('#consensus_taxonomy\tNCBI_tax_id\toriginal\n')
    empty, _ = load_profile(raw)
    assert empty.empty and list(empty) == ['taxid', 'N_original']


@pytest.mark.parametrize('body', [
    '#consensus_taxonomy\tsample\nEcoli\t1\n',
    '# call: motus profile -c\n#consensus_taxonomy\tNCBI_tax_id\tsample\n',
    '#consensus_taxonomy\tNCBI_tax_id\tsample\nEcoli [ref_mOTU_v3_1]\t562\t-0.1\n',
    '#consensus_taxonomy\tNCBI_tax_id\tsample\nEcoli [ref_mOTU_v3_1]\t562\tnan\n',
    '#consensus_taxonomy\tNCBI_tax_id\tsample\nEcoli [ref_mOTU_v3_1]\t562\t1.1\n',
    '#consensus_taxonomy\tNCBI_tax_id\tsample\nGenus\t562\t1\n',
    '#consensus_taxonomy\tNCBI_tax_id\tsample\na [ref_mOTU_v3_1]\t562\t.2\na [ref_mOTU_v3_1]\t562\t.2\n',
])
def test_motus3_rejects_incompatible_or_broken_profiles(tmp_path, body):
    raw = tmp_path/'bad.motus'; raw.write_text(body)
    with pytest.raises(ValueError):
        load_profile(raw, {'format': 'motus'})


@pytest.mark.parametrize('paired', [False, True])
def test_motus3_runner_configured_companions_and_private_paths(tmp_path, monkeypatch, paired):
    from samovar.profiles import run_motus
    import subprocess
    db = tmp_path/'database with spaces'; db.mkdir()
    (db/'db_mOTU_versions').write_text('motus\t3.1.0\n')
    (db/'db_mOTU_DB_NR.fasta.bwt').touch()
    tool = tmp_path/'configured tool'; tool.touch()
    monkeypatch.setattr('samovar.paths.resolve_executable', lambda *a, **k: str(tool))
    reads = tmp_path/'reads with spaces.fastq'
    reads.write_text('@r\nACGT\n+\nIIII\n')
    out = tmp_path/'output with spaces/raw.motus'
    calls = []
    def run(argv, *, cwd, env, check):
        calls.append(argv)
        assert check and argv[0] == str(tool)
        assert argv[argv.index('-t')+1] == '2'
        assert argv[argv.index('-g')+1] == '1'
        assert ('-f' in argv and '-r' in argv) == paired
        assert ('-s' in argv) != paired
        work = Path(cwd)
        assert (work/'db_mOTU').resolve() == db
        assert str(tool) in (work/'bin/bwa').read_text()
        assert env['PATH'].split(':')[0] == str(work/'bin')
        (work/'profile.tsv').write_text('#consensus_taxonomy\tNCBI_tax_id\ts\nEcoli [ref_mOTU_v3_1]\t562\t1\n')
    monkeypatch.setattr(subprocess, 'run', run)
    run_motus(reads, reads if paired else None, db, out, 's', 2, flags='-g 1')
    assert load_profile(out)[0].taxid.tolist() == ['562']
    assert not list((out.parent/'.tmp/motus').iterdir())
    reads.write_text(''); calls.clear()
    run_motus(reads, None, db, out, 's')
    assert load_profile(out)[0].empty and not calls
    with pytest.raises(SystemExit):
        run_motus(reads, None, db, out, 's', flags='-c')
    (db/'db_mOTU_DB_NR.fasta.bwt').unlink()
    with pytest.raises(FileNotFoundError):
        run_motus(reads, None, db, out, 's')


@pytest.mark.parametrize('outcome', ['failed', 'missing', 'broken'])
def test_motus3_native_errors_preserve_prior_output(tmp_path, monkeypatch, outcome):
    import subprocess
    from samovar.profiles import run_motus
    db = tmp_path/'db'; db.mkdir()
    (db/'db_mOTU_versions').write_text('motus\t3.1.0\n')
    (db/'db_mOTU_DB_NR.fasta.bwt').touch()
    tool = tmp_path/'motus'; tool.touch()
    monkeypatch.setattr('samovar.paths.resolve_executable', lambda *a, **k: str(tool))
    reads = tmp_path/'reads.fastq'; reads.write_text('@r\nACGT\n+\nIIII\n')
    out = tmp_path/'out'; out.write_text('previous result')
    def run(argv, **kwargs):
        if outcome == 'failed':
            raise subprocess.CalledProcessError(2, argv)
        if outcome == 'broken':
            (Path(kwargs['cwd'])/'profile.tsv').write_text('unexpected output')
    monkeypatch.setattr(subprocess, 'run', run)
    error = {'failed': subprocess.CalledProcessError, 'missing': FileNotFoundError, 'broken': ValueError}[outcome]
    with pytest.raises(error):
        run_motus(reads, None, db, out, 's')
    assert out.read_text() == 'previous result'
    assert not list((tmp_path/'.tmp/motus').iterdir())
    (db/'db_mOTU_versions').write_text('motus\t4.0.0\n')
    with pytest.raises(ValueError, match='version 3'):
        run_motus(reads, None, db, out, 's')
