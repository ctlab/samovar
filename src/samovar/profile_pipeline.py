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
    if stage == "initial":
        from samovar.exec_control import missing_profiler_reports

        missing = missing_profiler_reports(root)
        if missing:
            raise ValueError(
                "Missing profiler output, so the sample would be dropped from "
                "the combined profile: " + ", ".join(missing)
            )
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


def _drop_reprofiled_suffix(table, meta):
    """``A_reprofiled.csv`` is sample A, the same id as the initial reads."""
    rename = {
        col: col[: -len("_reprofiled")]
        for col in n_sample_columns(table)
        if str(col).endswith("_reprofiled")
    }
    if not rename:
        return table, meta
    table = table.rename(columns=rename)
    raw = dict((meta or {}).get("unassigned_by_sample") or {})
    if not raw:
        return table, meta
    meta = dict(meta)
    meta["unassigned_by_sample"] = {
        (key[: -len("_reprofiled")] if str(key).endswith("_reprofiled") else key): value
        for key, value in raw.items()
    }
    return table, meta


def stage_tables(root, stage):
    ann = stage_annotation(root, stage)
    tables = input_to_abundance_tables(ann)
    meta = dict(ann.profile_metadata)
    for name in tables:
        meta.setdefault(name, metadata({"rank": "taxid"}))
        if stage == "reprofiled":
            tables[name], meta[name] = _drop_reprofiled_suffix(tables[name], meta[name])
    return tables, meta


def _counts_table(columns):
    return (
        pd.DataFrame({
            name: pd.Series(dict(counts), dtype=float)
            for name, counts in columns.items()
        })
        .fillna(0)
        .rename_axis("taxid")
        .reset_index()
    )


def _read_folder(root, stage):
    """Reads that were sent to classifiers: trimmed folder when QC wrote one."""
    root = Path(root)
    trimmed = root / ("initial_trimmed" if stage == "initial" else "regenerated_trimmed")
    if trimmed.exists():
        return trimmed
    return root / ("initial" if stage == "initial" else "regenerated")


def empty_read_samples(root, stage):
    """FASTQ samples with no records. A missing per-read column is an empty profile."""
    from samovar.seqio import fastq_has_reads, find_fastq_mate, list_fastq_samples

    folder = _read_folder(root, stage)
    found = []
    for sample in list_fastq_samples(folder):
        path = find_fastq_mate(folder, sample)
        if path is None or not fastq_has_reads(path):
            found.append(sample)
    return found


def _fill_empty_read_samples(table, samples):
    missing = [f"N_{sample}" for sample in samples if f"N_{sample}" not in map(str, table.columns)]
    if not missing:
        return table
    table = table.copy()
    for col in missing:
        table[col] = 0.0
    return table


def truth_from_reads(root, stage):
    """Aggregate actual simulator truth, not the requested pre-simulation mass."""
    from collections import Counter
    from samovar.parse_annotators import extract_true_taxid
    from samovar.seqio import list_fastq_samples, find_fastq_mate, iter_fastq_records
    root = Path(root)
    folder = _read_folder(root, stage)
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
    table = _counts_table(columns)
    return table, metadata({
        "rank": "taxid",
        "normalization": "all",
        "source": "genome-based read headers",
        "ground_truth": "reads",
    })


def _truth_from_true_column(root, stage):
    from samovar.abundance import annotation_to_abundance_tables

    frame = stage_annotation(root, stage).DataFrame
    if frame.empty or "true" not in frame or not frame["true"].notna().all():
        return None
    labels = frame["true"].astype(str)
    if not (labels.str.fullmatch(r"[0-9]+").all() and labels.ne("0").any()):
        return None
    truth = frame[[c for c in ("sample", "true") if c in frame]].rename(columns={"true": "taxID_truth"})
    table = annotation_to_abundance_tables(truth)["truth"]
    return table, metadata({
        "rank": "taxid",
        "normalization": "all",
        "source": "per-read true column",
        "ground_truth": "per-read",
    })


def _read_folder(root, stage):
    root = Path(root)
    folder = root / ("initial_trimmed" if stage == "initial" else "regenerated_trimmed")
    if not folder.is_dir():
        folder = root / ("initial" if stage == "initial" else "regenerated")
    return folder


def _truth_lookup(pairs):
    """Same read-id matching as the annotation combiner."""
    exact = {}
    prefixes = {}
    for key, tax in pairs:
        if not key or not tax or tax == "0":
            continue
        exact[key] = tax
        if len(key) > 2 and key[-2] == "/":
            exact[key[:-2]] = tax
        if "/" not in key and len(key) >= 4:
            prefixes[key] = tax

    def lookup(seq):
        if seq in exact:
            return exact[seq]
        token = seq.split()[0]
        if token in exact:
            return exact[token]
        if len(token) > 2 and token[-2] == "/" and token[:-2] in exact:
            return exact[token[:-2]]
        best = ""
        for key in prefixes:
            if len(key) <= len(best) or len(key) > len(token) or not token.startswith(key):
                continue
            best = key
        return prefixes.get(best, "")

    return lookup


def profile_from_truth_table(root, table_path, stage="initial"):
    """Count a seq→taxid table onto FASTQ samples. None when no FASTQ is present."""
    from collections import Counter

    from samovar.ground_truth import iter_truth_pairs
    from samovar.seqio import fastq_record_id, find_fastq_mate, iter_fastq_records, list_fastq_samples

    path = Path(table_path)
    pairs = [(seq, tax) for seq, tax in iter_truth_pairs(path) if tax and tax != "0"]
    if not pairs:
        raise ValueError(f"{path} is not a seq→taxid ground-truth table")
    lookup = _truth_lookup(pairs)
    folder = _read_folder(root, stage)
    samples = list_fastq_samples(folder)
    if not samples:
        return None
    columns = {}
    matched = 0
    for sample in samples:
        counts = Counter()
        for header, _seq, _plus, _qual in iter_fastq_records(find_fastq_mate(folder, sample)):
            tax = lookup(fastq_record_id(header))
            if not tax:
                continue
            counts[tax] += 1
            matched += 1
        columns["N_" + sample] = counts
    if matched == 0:
        raise ValueError(f"No read ids under {folder} match the seq→taxid table {path}")
    return _counts_table(columns), metadata({
        "rank": "taxid",
        "normalization": "all",
        "source": "per-read seq→taxid table",
        "ground_truth": "per-read-table",
    })


def ground_truths(root, stage, config):
    """Every supplied ground truth, with genome-based reads first."""
    read_stage = "initial" if stage == "reprofiled" else stage
    if stage == "regenerated":
        aggregated = _truth_from_true_column(root, "regenerated")
        if aggregated:
            return [("per-read", *aggregated)]
        headers = truth_from_reads(root, "regenerated")
        if headers:
            table, meta = headers
            return [("reads", table, {**meta, "ground_truth": "reads"})]
        return []
    found = []
    headers = truth_from_reads(root, read_stage)
    if headers:
        table, meta = headers
        found.append(("reads", table, {**meta, "ground_truth": "reads"}))
    table_path = str(
        config.get("initial_ground_truth_table") or config.get("truth_table") or ""
    ).strip()
    if table_path and read_stage == "initial":
        mapped = profile_from_truth_table(root, table_path, "initial")
        if mapped:
            found.append(("per-read-table", *mapped))
    if stage in {"initial", "reprofiled"} and config.get("profile_truth"):
        desc = config["profile_truth"]
        table, meta = load_profile(desc["path"], desc)
        found.append(("profile_truth", table, {
            **meta,
            "ground_truth": "profile_truth",
            "source": meta.get("source") or "profile_truth",
        }))
    if not found:
        aggregated = _truth_from_true_column(root, read_stage)
        if aggregated:
            found.append(("per-read", *aggregated))
    return found


def ground_truth(root, stage, config):
    """Primary ground truth. Genome-based reads come before supplied alternatives."""
    found = ground_truths(root, stage, config)
    if not found:
        return None
    return found[0][1], found[0][2]


def _taxdump_present(taxdump) -> bool:
    from samovar.kreport import resolve_taxdump_files

    try:
        resolve_taxdump_files(taxdump or None)
    except FileNotFoundError:
        return False
    return True


def align_truth(pred, pm, truth, tm, root, config):
    if pm["taxonomy"] != tm["taxonomy"]:
        raise ValueError(f"Incompatible profile taxonomy: {pm['taxonomy']} vs {tm['taxonomy']}; supply an explicit taxonomy mapping")
    taxdump = config.get("taxdump") or None
    # Per-read tables and read-header truth are both rank taxid. Species
    # rollup needs a taxdump; without one, compare the taxids as stored.
    if pm["rank"] == "taxid" and tm["rank"] == "taxid" and not _taxdump_present(taxdump):
        rank = "taxid"
    else:
        rank = pm["rank"] if pm["rank"] != "taxid" else "species"
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
    native, native_meta = read_profiles(root / f"{stage}_profiles")
    tables, metas = stage_tables(root, stage)
    read_stage = "initial" if stage == "reprofiled" else stage
    empty = empty_read_samples(root, read_stage)
    for name in tables:
        tables[name] = _fill_empty_read_samples(tables[name], empty)
    for name, table in native.items():
        if name not in tables:
            tables[name] = table
            metas[name] = native_meta[name]
    gts = ground_truths(root, stage, config)
    if not tables or (not gts and not (native or config.get("profile_truth"))):
        return False
    from samovar.scores import profile_scores
    rows = []
    if not gts:
        for name in tables:
            rows.append(pd.DataFrame([dict(annotator=name, reason="no profile gold standard available")]))
    for gt_name, truth, tm in gts:
        for name, table in tables.items():
            pred, pm, aligned, aligned_meta = align_truth(table, metas[name], truth, tm, root, config)
            result = profile_scores(pred, aligned, pm, aligned_meta)
            result.insert(0, "ground_truth", gt_name)
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
    opal_gt = next(((truth, tm) for _name, truth, tm in gts if tm.get("kind") == "abundance"), None)
    if config.get("viz", {}).get("opal") and opal_gt is not None:
        from samovar.cami_profile import write_cami_profile
        from samovar.opal import run_opal
        profile_files = []
        gold_files = []
        for name, table in tables.items():
            pred, pm, truth, tm = align_truth(table, metas[name], *opal_gt, root, config)
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
    model_name = "profile_joint.json" if method == "profile_joint" else "profile_ratio.json"
    cfg = {**config, "reference": reference, "reference_truth": reference_truth,
           "run_root": str(root), "profile_model": str(root / "reprofiled_profiles" / model_name)}
    return run_export(method, stage_annotation(root, "initial"), root / "reprofiled_profiles", cfg)
