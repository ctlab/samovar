"""Profile boundaries of existing combine, abundance, scoring and export stages.

Stage scheduling/checkpoints remain in PipelineConfig and Snakemake. These
functions only load/materialize the independent sample-level representation.
"""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from samovar.abundance import input_to_abundance_tables, n_sample_columns
from samovar.profiles import (at_rank, change_basis, load_profile, metadata, read_counts,
                              read_profiles, reference_lengths, write_profiles)


def collect(root, stage, config):
    root = Path(root)
    grouped, metas = {}, {}
    sources = [(p, {}) for p in sorted((root / f"{stage}_reports").glob("*.profile.csv"))]
    if stage == "initial":
        sources += [(Path(info["path"]), {**info, "identity": name})
                    for name, info in (config.get("profile_inputs") or {}).items()]
    for path, options in sources:
        table, meta = load_profile(path, options)
        name = meta.get("identity")
        if not name:
            raise ValueError(f"Profile lacks an explicit identity: {path}")
        if name in grouped:
            previous = metas[name]
            for key in ("unit", "basis", "kind", "taxonomy", "rank", "normalization", "database"):
                if previous.get(key) != meta.get(key):
                    raise ValueError(f"Inconsistent {key} for profile {name}")
            overlap = set(n_sample_columns(table)) & set(n_sample_columns(grouped[name]))
            if overlap:
                raise ValueError(f"Duplicate samples for profile {name}: {overlap}")
            grouped[name] = grouped[name].merge(table, on="taxid", how="outer").fillna(0)
        else:
            grouped[name], metas[name] = table, meta
        metas[name].setdefault("unassigned_by_sample", {}).update(meta.get("unassigned_by_sample", {}))
    if grouped:
        write_profiles(root / f"{stage}_profiles", grouped, metas)
    return grouped, metas


def stage_annotation(root, stage):
    from samovar.parse_annotators import Annotation
    root = Path(root)
    obj = Annotation.from_annotation_dir(root / f"{stage}_annotations")
    return obj


def stage_tables(root, stage):
    ann = stage_annotation(root, stage)
    tables = input_to_abundance_tables(ann)
    meta = dict(ann.profile_metadata)
    for name in tables:
        meta.setdefault(name, metadata({"rank": "taxid"}))
    return tables, meta


def truth_from_reads(root, stage):
    """Aggregate actual simulator truth, not the requested pre-simulation mass."""
    from collections import Counter
    from samovar.parse_annotators import extract_true_taxid
    from samovar.seqio import list_fastq_samples, find_fastq_mate, iter_fastq_records
    root = Path(root)
    # Truth describes the reads actually sent to classifiers (after QC).
    folder = root / ("initial_trimmed" if stage == "initial" else "regenerated_trimmed")
    if not folder.exists():
        folder = root / ("initial" if stage == "initial" else "regenerated")
    columns = {}
    for sample in list_fastq_samples(folder):
        path = find_fastq_mate(folder, sample)
        counts = Counter()
        for rec in iter_fastq_records(path):
            tax = extract_true_taxid(rec[0])
            if not tax or tax == "0":
                if stage == "initial":
                    return None  # real reads do not imply an empty gold standard
                raise ValueError(f"Regenerated read has missing taxid truth: {rec[0].strip()}")
            counts[tax] += 1
        columns["N_" + sample] = counts
    if not columns:
        return None
    table = pd.DataFrame({name: pd.Series(dict(counts), dtype=float) for name, counts in columns.items()}).fillna(0).rename_axis("taxid").reset_index()
    return table, metadata({"rank": "taxid", "normalization": "all", "source": "actual simulator read headers"})


def ground_truth(root, stage, config):
    if stage in {"initial", "reprofiled"} and config.get("profile_truth"):
        desc = config["profile_truth"]
        return load_profile(desc["path"], desc)
    from samovar.abundance import annotation_to_abundance_tables
    read_stage = "initial" if stage == "reprofiled" else stage
    frame = stage_annotation(root, read_stage).DataFrame
    if not frame.empty and "true" in frame and frame["true"].notna().all():
        labels = frame["true"].astype(str)
        if labels.str.fullmatch(r"[0-9]+").all() and labels.ne("0").any():
            truth = frame[[c for c in ("sample", "true") if c in frame]].rename(columns={"true": "taxID_truth"})
            table = annotation_to_abundance_tables(truth)["truth"]
            return table, metadata({"rank": "taxid", "normalization": "all", "source": "per-read true column"})
    return truth_from_reads(root, "initial" if stage == "reprofiled" else stage)


def align_truth(pred, pm, truth, tm, root, config):
    if pm["taxonomy"] != tm["taxonomy"]:
        raise ValueError(f"Incompatible profile taxonomy: {pm['taxonomy']} vs {tm['taxonomy']}; supply an explicit taxonomy mapping")
    rank = pm["rank"] if pm["rank"] != "taxid" else "species"
    taxdump = config.get("taxdump")
    pred, pm = at_rank(pred, pm, rank, taxdump)
    truth, tm = at_rank(truth, tm, rank, taxdump)
    if tm["kind"] != "presence" and pm["normalization"] == "classified" and tm["normalization"] == "all":
        from samovar.profiles import fractions
        truth = fractions(truth, tm).query("taxid != '0'").copy()
        cols = n_sample_columns(truth)
        truth[cols] = truth[cols].div(truth[cols].sum().replace(0, float("nan"))).fillna(0)
        tm = {**tm, "unit": "fraction", "normalization": "classified", "normalization_conversion": "all to classified"}
    if tm["kind"] != "presence" and tm["basis"] != pm["basis"]:
        tids = set(truth.loc[truth[n_sample_columns(truth)].sum(axis=1) > 0, "taxid"]) - {"0"}
        lengths = reference_lengths(tids, Path(root) / "genomes", config)
        truth, tm = change_basis(truth[truth.taxid != "0"], tm, pm["basis"], lengths)
    return pred, pm, truth, tm


def score(root, stage, config):
    """Write separate profile metrics and MultiQC, leaving per-read scores alone."""
    root = Path(root)
    config = {**config.get("profile_settings", {}), **config}
    native, _ = read_profiles(root / f"{stage}_profiles")
    if not native and not config.get("profile_truth"):
        return False
    from samovar.scores import profile_scores
    tables, metas = stage_tables(root, stage)
    gt = ground_truth(root, stage, config)
    rows = []
    for name, table in tables.items():
        if gt is None:
            rows.append(pd.DataFrame([dict(annotator=name, reason="no profile gold standard available")]))
            continue
        pred, pm, truth, tm = align_truth(table, metas[name], *gt, root, config)
        result = profile_scores(pred, truth, pm, tm)
        result.insert(0, "annotator", name)
        rows.append(result)
    dest = root / f"{stage}_annotations_plots"
    dest.mkdir(parents=True, exist_ok=True)
    result = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    result.to_csv(dest / "profile_quality_scores.csv", index=False)
    data = {str(i): {k: (None if pd.isna(v) else v) for k, v in row.items()}
            for i, row in result.iterrows()}
    content = dict(id=f"{stage}_profile_scores", section_name="Taxonomic profiles",
                   description="Sample-level metrics; no per-read assignments inferred.", plot_type="table", data=data)
    (dest / "profile_quality_scores_mqc.json").write_text(json.dumps(content, indent=2, allow_nan=False) + "\n")
    if config.get("viz", {}).get("opal") and gt is not None and gt[1]["kind"] == "abundance":
        from samovar.cami_profile import write_cami_profile
        from samovar.opal import run_opal
        profile_files = []
        gold_files = []
        for name, table in tables.items():
            pred, pm, truth, tm = align_truth(table, metas[name], *gt, root, config)
            from samovar.profiles import fractions
            pred, truth = fractions(pred, pm), fractions(truth, tm)
            for col in n_sample_columns(pred):
                folder = dest / "profile_opal" / col[2:]
                gold = write_cami_profile(truth.set_index("taxid")[col].to_dict(), folder / "gold.profile", col[2:], rank=pm["rank"], total=1, taxonomy_id=pm["taxonomy"])
                path = write_cami_profile(pred.set_index("taxid")[col].to_dict(), folder / f"{name}.profile", col[2:], rank=pm["rank"], total=1, taxonomy_id=pm["taxonomy"])
                if gold and path:
                    profile_files.append((name, path))
                    gold_files.append(gold)
        # Each sample has an independent gold/profile pair and benchmark directory.
        for (name, path), gold in zip(profile_files, gold_files):
            run_opal(gold, [(name, path)], path.parent / name, rank="species")
    return True


def materialize(root, config):
    """Compile raw typed profiles to read counts before existing table regeneration."""
    from samovar.abundance import collect_observed_abundance, write_abundance_dir
    root = Path(root)
    raw, metas = read_profiles(root / "initial_profiles")
    tables = collect_observed_abundance(root)
    provenance = {}
    for name, table in raw.items():
        info = metas[name]
        if info["kind"] == "presence":
            raise ValueError("Presence input cannot define a quantitative regeneration community")
        if info["taxonomy"] != "ncbi":
            raise ValueError("Profile regeneration needs an explicit mapping to NCBI reference taxids")
        table, info = at_rank(table, info, "species", config.get("taxdump"))
        table = table[table.taxid != "0"].copy()
        tids = set(table.loc[table[n_sample_columns(table)].sum(axis=1) > 0, "taxid"])
        lengths = reference_lengths(tids, root / "genomes", config) if info["basis"] != "reads" else {}
        if info["unit"] == "counts" and info["basis"] == "reads":
            converted, cm = table, info
        else:
            converted, cm = read_counts(table, info, int(config.get("N_reads", 1000)), lengths)
        tables[name] = converted
        provenance[name] = cm
    write_abundance_dir(root / "initial_abundance", tables)
    if provenance:
        (root / "initial_abundance" / ".profile_conversion.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return tables


def correct(root, config):
    from samovar.abundance_correctors import run_export
    root = Path(root)
    config = {**config.get("profile_settings", {}), **config}
    initial, _ = read_profiles(root / "initial_profiles")
    if not initial:
        return
    method = config.get("profile_corrector", "profile_ratio")
    if method == "none":
        tables, metas = read_profiles(root / "initial_profiles")
        return write_profiles(root / "reprofiled_profiles", tables, metas)
    reference = stage_annotation(root, "regenerated")
    reference_truth = ground_truth(root, "regenerated", config)
    if reference_truth is None:
        raise ValueError("Profile correction needs actual regenerated truth")
    cfg = {**config, "reference": reference, "reference_truth": reference_truth,
           "run_root": str(root), "profile_model": str(root / "reprofiled_profiles" / "profile_ratio.json")}
    return run_export(method, stage_annotation(root, "initial"), root / "reprofiled_profiles", cfg)
