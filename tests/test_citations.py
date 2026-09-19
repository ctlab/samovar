import json
from pathlib import Path

from samovar.citations import (
    DOI_RECORDS,
    STATIC_RECORDS,
    cite_dir,
    mapping_for_tools,
    missing_tools,
    rebuild,
)
from samovar.main_config import TOOL_GROUP_BY_NAME


def test_every_builtin_tool_is_in_the_map():
    mapping = mapping_for_tools()
    for name in TOOL_GROUP_BY_NAME:
        assert name in mapping
    assert "samovar" in mapping
    assert "sparsedossa2" in mapping


def test_listed_bib_files_are_known_records():
    known = set(DOI_RECORDS) | set(STATIC_RECORDS)
    for name, files in mapping_for_tools().items():
        for filename in files:
            assert filename in known, f"{name} -> {filename}"


def test_aliases_share_canonical_papers():
    mapping = mapping_for_tools()
    assert mapping["metaphlan4"] == mapping["metaphlan"]
    assert mapping["art_illumina"] == mapping["art"]
    assert mapping["opal.py"] == mapping["opal"]
    assert mapping["DAS_Tool"] == mapping["dastool"]
    assert mapping["nanosim3"] == mapping["nanosim"]


def test_rebuild_offline_writes_json(tmp_path):
    mapping = rebuild(tmp_path, online=False, harvest=False)
    dest = cite_dir(tmp_path)
    assert (dest / "citations.json").is_file()
    assert (dest / "missing.json").is_file()
    disk = json.loads((dest / "citations.json").read_text())
    assert disk["kraken2"] == mapping["kraken2"]
    assert (dest / "kraken2_wood2019.bib").is_file()
    assert (dest / "sklearn_pedregosa2011.bib").is_file()
    assert (dest / "samovar_chechenina2023.bib").is_file()
    for filename in mapping["kraken2"]:
        assert (dest / filename).is_file()
    missing = json.loads((dest / "missing.json").read_text())
    assert missing == missing_tools(mapping)
    for name in missing:
        assert disk[name] == []


def test_missing_tools_are_runtimes_and_seqtk():
    missing = set(missing_tools())
    assert "seqtk" in missing
    assert "bash" in missing
    assert "python" in missing
    assert "g++" in missing
    assert "kraken2" not in missing
    assert "iss" not in missing
    assert "camisim" not in missing


def test_repo_cite_tree_matches_bundled_map():
    root = Path(__file__).resolve().parents[1]
    dest = cite_dir(root)
    disk = json.loads((dest / "citations.json").read_text())
    expected = mapping_for_tools()
    bundled = {k: [f for f in v if not str(f).endswith("_harvest.bib")] for k, v in disk.items()}
    assert bundled == expected
    for files in expected.values():
        for filename in files:
            assert (dest / filename).is_file()
            text = (dest / filename).read_text()
            assert "@" in text


def test_write_used_citations_from_prepare_config(tmp_path):
    import pytest

    from samovar.citations import selected_prepare_tools, write_used_citations
    from samovar.config import AnnotatorConfig, PipelineConfig

    cfg = PipelineConfig()
    cfg.output_dir = str(tmp_path)
    cfg.annotators = [
        AnnotatorConfig(run_name="k2", type="kraken2", db_path="/db", cmd="kraken2"),
        AnnotatorConfig(run_name="kj", type="kaiju", db_path="/db", cmd="kaiju"),
    ]
    cfg.reads_generator = "iss"
    cfg.regeneration_modes = ["direct"]
    cfg.regeneration_mode = "direct"
    cfg.run_multiqc = False
    cfg.export_corrector = "off"
    cfg.scoring_tools = []
    names = selected_prepare_tools(cfg)
    assert "kraken2" in names
    assert "kaiju" in names
    assert "iss" in names
    assert "snakemake" in names
    cfg.annotators.append(
        AnnotatorConfig(
            run_name="asm",
            type="assembly",
            db_path=".",
            cmd="assembly",
            extra="--assembler megahit --gene-caller prodigal --binner metabat2 --aligner minimap2",
        )
    )
    names = selected_prepare_tools(cfg)
    assert "megahit" in names
    assert "prodigal" in names
    assert "metabat2" in names
    assert "minimap2" in names
    with pytest.warns(UserWarning, match="No BibTeX citation"):
        path = write_used_citations(cfg, tmp_path)
    text = path.read_text(encoding="utf-8")
    assert (tmp_path / "used_citations.bib").is_file()
    assert (tmp_path / ".log" / "configs" / "used_citations.bib").is_file()
    assert "Kraken 2" in text or "Wood_" in text
    assert "Kaiju" in text or "Menzel" in text
    assert "InSilicoSeq" in text or "Gourl" in text
    assert "Snakemake" in text or "snakemake" in text.lower() or "Koster" in text
    with pytest.warns(UserWarning, match="No BibTeX citation"):
        again = write_used_citations(cfg, tmp_path).read_text(encoding="utf-8")
    assert again == text

