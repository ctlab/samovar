"""Sample-level taxonomic profiles: typed abundance tables, never invented reads.

The numeric representation is the existing taxid x N_<sample> table. Metadata
travels beside it so organism abundance cannot silently become read counts.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from samovar.abundance import n_sample_columns

FORMAT_VERSION = 1
MANIFEST = ".profiles.json"
DEFAULT_METADATA = dict(kind="abundance", unit="counts", basis="reads",
                        taxonomy="ncbi", rank="species", normalization="classified")


def metadata(settings=None):
    out = {**DEFAULT_METADATA, **dict(settings or {})}
    choices = dict(kind={"abundance", "presence"}, unit={"counts", "fraction", "percent", "binary"},
                   basis={"reads", "organisms", "presence"}, taxonomy={"ncbi", "gtdb"},
                   normalization={"classified", "all"})
    for key, allowed in choices.items():
        if out[key] not in allowed:
            raise ValueError(f"Invalid profile {key}: {out[key]!r}")
    if out["kind"] == "presence":
        out.update(unit="binary", basis="presence")
    elif out["unit"] == "binary" or out["basis"] == "presence":
        raise ValueError("Binary profiles require kind=presence")
    if not isinstance(out["rank"], str) or not out["rank"]:
        raise ValueError("Profile rank is required")
    return out


def validate(table, settings=None):
    """Strict boundary validation; do not coerce broken cells into absence."""
    meta = metadata(settings)
    work = table.copy()
    if "taxid" not in work or not n_sample_columns(work):
        raise ValueError("Profile requires taxid and explicit N_<sample> columns, including empty profiles")
    if work.columns.duplicated().any() or work["taxid"].isna().any():
        raise ValueError("Duplicate profile columns or missing taxid")
    work["taxid"] = work["taxid"].astype(str).str.strip()
    if work["taxid"].eq("").any() or work["taxid"].duplicated().any():
        raise ValueError("Empty or duplicate profile taxid")
    if meta["taxonomy"] == "ncbi" and not work["taxid"].map(lambda s: s.isdigit() or s.startswith("unmapped:")).all():
        raise ValueError("NCBI profiles need integer taxids; preserve unresolvable native taxa as unmapped:<native-id>")
    cols = n_sample_columns(work)
    values = work[cols].apply(pd.to_numeric, errors="raise").astype(float)
    if not np.isfinite(values.to_numpy()).all() or (values < 0).any().any():
        raise ValueError("Profile values must be finite and nonnegative")
    if meta["kind"] == "presence" and not values.isin([0, 1]).all().all():
        raise ValueError("Presence gold standard accepts only 0 or 1")
    if meta["unit"] in {"fraction", "percent"}:
        limit = 1 if meta["unit"] == "fraction" else 100
        if (values.sum() > limit + 1e-5).any():
            raise ValueError(f"Profile {meta['unit']} column mass exceeds {limit}; check ranks/units")
    work = work[["taxid", *cols]]
    work[cols] = values
    return work, meta


def load_profile(path, settings=None):
    """Read canonical CSV/TSV, CAMI or native MetaPhlAn at one declared rank."""
    path = Path(path)
    sidecar = Path(str(path) + ".profile.json")
    saved = json.loads(sidecar.read_text()) if sidecar.is_file() else {}
    meta = metadata({**saved, **dict(settings or {})})
    fmt = meta.get("format", "auto")
    with path.open() as handle:
        lines = handle.readlines()
    if fmt == "auto":
        fmt = "cami" if any(x.startswith("@@TAXID") for x in lines) else (
            "metaphlan" if any("clade_name" in x and "relative_abundance" in x for x in lines) else "abundance")
    if fmt == "abundance":
        table = pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",", dtype={"taxid": str})
    elif fmt in {"metaphlan", "cami"}:
        rows = []
        sample = meta.get("sample", path.stem)
        declared_samples = []
        rank_codes = dict(superkingdom="k", kingdom="k", phylum="p", **{"class": "c"},
                          order="o", family="f", genus="g", species="s")
        for line in lines:
            if line.startswith("@SampleID:") and "sample" not in meta:
                sample = line.split(":", 1)[1].strip()
                declared_samples.append(str(sample))
            if line.startswith("@TaxonomyID:"):
                system = line.split(":", 1)[1].strip().lower()
                if settings and settings.get("taxonomy", system) != system:
                    raise ValueError("CAMI taxonomy disagrees with profile metadata")
                meta["taxonomy"] = system
            if not line.strip() or line.startswith(("#", "@")):
                continue
            fields = line.rstrip().split("\t")
            if fmt == "cami":
                if len(fields) < 4:
                    raise ValueError(f"Malformed CAMI profile row: {line!r}")
                if fields[1] == meta["rank"]:
                    rows.append((str(sample), fields[0], float(fields[-1])))
            else:
                if len(fields) < 3:
                    raise ValueError(f"Malformed MetaPhlAn profile row: {line!r}")
                leaf = fields[0].split("|")[-1]
                if leaf == "UNCLASSIFIED":
                    meta.setdefault("unassigned_by_sample", {})[str(sample)] = float(fields[2]) / 100
                    meta["normalization"] = "all"
                elif leaf.startswith(rank_codes.get(meta["rank"], "?") + "__"):
                    tax = fields[1].split("|")[-1]
                    if not tax.isdigit() or int(tax) <= 0:
                        tax = "unmapped:" + fields[0]
                    rows.append((str(sample), tax, float(fields[2])))
        samples = list(dict.fromkeys([*declared_samples, *(row[0] for row in rows)])) or [str(sample)]
        table = pd.DataFrame(rows, columns=["sample", "taxid", "value"])
        if table.empty:
            table = pd.DataFrame(columns=["taxid", *("N_" + s for s in samples)])
        else:
            table = table.pivot_table(index="taxid", columns="sample", values="value", aggfunc="sum", fill_value=0)
            table = table.reindex(columns=samples, fill_value=0)
            table.columns = ["N_" + str(s) for s in table.columns]
            table = table.reset_index()
        meta.update(unit="percent", basis="organisms" if fmt == "metaphlan" else meta["basis"])
    else:
        raise ValueError(f"Unknown profile format: {fmt}")
    meta.pop("format", None)
    return validate(table, meta)


def write_profiles(dest, tables, meta):
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    saved = {}
    for name, table in tables.items():
        if not name or Path(name).name != name or name in {".", ".."}:
            raise ValueError(f"Invalid profile identity: {name!r}")
        table, info = validate(table, meta.get(name))
        table.to_csv(dest / f"{name}.csv", index=False)
        saved[name] = info
    (dest / MANIFEST).write_text(json.dumps(dict(version=FORMAT_VERSION, profiles=saved), indent=2) + "\n")
    return dest


def read_profiles(source):
    source = Path(source)
    manifest = source / MANIFEST
    if not manifest.is_file():
        return {}, {}
    saved = json.loads(manifest.read_text())
    if saved.get("version") != FORMAT_VERSION or not isinstance(saved.get("profiles"), dict):
        raise ValueError(f"Unsupported/corrupt profile manifest: {manifest}")
    tables, metas = {}, {}
    for name, info in saved["profiles"].items():
        if Path(name).name != name:
            raise ValueError("Invalid profile manifest identity")
        tables[name], metas[name] = load_profile(source / f"{name}.csv", info)
    return tables, metas


def fractions(table, info):
    table, info = validate(table, info)
    if info["kind"] != "abundance":
        raise ValueError("Presence has no quantitative abundance")
    cols = n_sample_columns(table)
    if info["normalization"] == "classified":
        table = table[table.taxid != "0"].copy()
    values = table[cols]
    # All-read fractions may legitimately sum to less than one.
    if info["unit"] == "percent":
        table[cols] = values / 100
    elif info["unit"] == "counts":
        table[cols] = values.div(values.sum().replace(0, np.nan)).fillna(0)
    return table


def change_basis(table, info, basis, lengths):
    table = fractions(table, info)
    result = {**metadata(info), "unit": "fraction"}
    if result["basis"] != basis:
        if result["normalization"] != "classified":
            mass = table[n_sample_columns(table)].sum()
            unassigned = table.loc[table.taxid == "0", n_sample_columns(table)].sum()
            if ((mass > 0) & ~np.isclose(mass, 1)).any() or (unassigned > 0).any():
                raise ValueError("Basis conversion cannot resolve genome length for unassigned mass; use a classified-normalized profile")
            table = table[table.taxid != "0"].copy()
        missing = set(table.loc[table[n_sample_columns(table)].sum(axis=1) > 0, "taxid"]) - set(lengths)
        if missing:
            raise ValueError(f"Missing selected reference lengths for profile taxa: {sorted(missing)}")
        weights = table["taxid"].map(lengths).fillna(1).astype(float)
        if (weights <= 0).any():
            raise ValueError("Selected reference genome length must be positive")
        cols = n_sample_columns(table)
        table[cols] = table[cols].mul(weights if basis == "reads" else 1 / weights, axis=0)
        table[cols] = table[cols].div(table[cols].sum().replace(0, np.nan)).fillna(0)
        result["basis"] = basis
        result["reference_lengths"] = dict(lengths)
        result["approximation"] = "selected-reference genome-length weighting"
    return table, result


def reference_lengths(taxids, genome_dir, config=None):
    """Use the same genome resolver/selection as read regeneration."""
    from samovar.table2iss import _resolve_genomes_for_taxids
    from samovar.seqio import fasta_handle
    cfg = dict(config or {})
    genome_dir = Path(genome_dir)
    genome_dir.mkdir(parents=True, exist_ok=True)
    paths = _resolve_genomes_for_taxids(taxids, str(genome_dir), cfg.get("email", ""), True,
                                      gzip_genomes=cfg.get("gzip_genomes", True),
                                      reannotation_level=cfg.get("reannotation_level", "taxid"),
                                      max_genome_mb=cfg.get("max_genome_mb"),
                                      genome_skip_list=cfg.get("genome_skip_list"))
    sizes = {}
    manifest = genome_dir / ".profile_references.json"
    saved = json.loads(manifest.read_text()) if manifest.is_file() else {}
    for tax, path in paths.items():
        with fasta_handle(path) as handle:
            sizes[tax] = sum(len(line.strip()) for line in handle if not line.startswith(">"))
        stat = Path(path).stat()
        saved[tax] = dict(path=str(Path(path).resolve()), length=sizes[tax],
                          size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    missing = set(taxids) - set(sizes)
    if missing:
        raise ValueError(f"No reference genomes for profile taxa: {sorted(missing)}")
    manifest.write_text(json.dumps(saved, indent=2) + "\n")
    return sizes


def read_counts(table, info, depth, lengths=None):
    """Convert a typed profile to simulator fragment counts with exact totals."""
    from samovar.table2iss import _scale_amounts
    if int(depth) < 0:
        raise ValueError("Profile generation depth must be nonnegative")
    table, converted = change_basis(table, info, "reads", lengths or {})
    for col in n_sample_columns(table):
        values = table[col].to_numpy(float)
        table[col] = _scale_amounts(values.tolist(), int(depth), fractional=True) if values.sum() > 0 else 0
    return table, {**converted, "unit": "counts", "depth": int(depth)}


def validate_profile_settings(config):
    inputs = config.get("profile_inputs") or {}
    if not isinstance(inputs, dict):
        raise ValueError("profile_inputs must map profile identities to descriptors")
    def descriptor(value):
        if not isinstance(value, dict) or not value.get("path"):
            raise ValueError("Profile descriptor requires a path and explicit metadata")
        result = dict(value)
        result["path"] = str(Path(result["path"]).expanduser().resolve())
        for key in ("kind", "unit", "basis", "taxonomy", "rank", "normalization"):
            if key not in result:
                raise ValueError(f"Profile descriptor requires explicit {key}")
        metadata(result)
        if not Path(result["path"]).is_file():
            raise FileNotFoundError(result["path"])
        load_profile(result["path"], result)
        return result
    inputs = {str(name): descriptor(value) for name, value in inputs.items()}
    truth = descriptor(config["profile_truth"]) if config.get("profile_truth") else None
    corrector = config.get("profile_corrector", "profile_ratio")
    if corrector not in {"profile_ratio", "none"}:
        from samovar.abundance_correctors import require_known_export
        require_known_export(corrector)
    return inputs, truth, corrector


def at_rank(table, info, rank, taxdump=None):
    """Project descendants upward using the existing taxonomy implementation."""
    table, info = validate(table, info)
    if info["rank"] == rank:
        return table, info
    from samovar.kreport import KReportTaxonomy
    tree = KReportTaxonomy.from_taxdump(taxdump, taxonomy=info["taxonomy"])
    mapped = []
    for token in table["taxid"]:
        node = tree.resolve(token)
        if node == 0 and token != "0":
            raise ValueError(f"Unresolved {info['taxonomy']} profile taxon {token}")
        mapped.append(str(next((p for p in tree.ancestors_inclusive(node) if tree.rank_of(p) == rank), 0)))
    table["taxid"] = mapped
    table = table.groupby("taxid", as_index=False)[n_sample_columns(table)].sum()
    return table, {**info, "rank": rank}


def compatible(pred, truth):
    for field in ("taxonomy", "rank"):
        if pred[field] != truth[field]:
            raise ValueError(f"Incompatible profile {field}: {pred[field]} vs {truth[field]}")
    if truth["kind"] != "presence":
        for field in ("basis", "normalization"):
            if pred[field] != truth[field]:
                raise ValueError(f"Incompatible profile {field}: {pred[field]} vs {truth[field]}")


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["ingest", "collect", "materialize", "correct"])
    parser.add_argument("--source")
    parser.add_argument("--output")
    parser.add_argument("--sample")
    parser.add_argument("--identity")
    parser.add_argument("--database", default="")
    parser.add_argument("--format", default="auto")
    parser.add_argument("--output-dir")
    parser.add_argument("--stage", default="initial")
    parser.add_argument("--config")
    args = parser.parse_args(argv)
    if args.action == "ingest":
        if args.format == "abundance" and not Path(str(args.source) + ".profile.json").is_file():
            raise ValueError("Profiler contract requires an output .profile.json sidecar")
        table, meta = load_profile(args.source, {"format": args.format, "sample": args.sample})
        cols = n_sample_columns(table)
        if len(cols) != 1:
            raise ValueError("Profiler execution must produce one sample")
        table = table.rename(columns={cols[0]: "N_" + args.sample})
        meta.update(identity=args.identity, database=args.database, source=str(args.source))
        table.to_csv(args.output, index=False)
        Path(args.output + ".profile.json").write_text(json.dumps(meta, indent=2) + "\n")
        return 0
    import yaml
    from samovar import profile_pipeline
    cfg = yaml.safe_load(Path(args.config).read_text()) if args.config else {}
    if args.action == "collect":
        profile_pipeline.collect(args.output_dir, args.stage, cfg or {})
    elif args.action == "materialize":
        profile_pipeline.materialize(args.output_dir, cfg or {})
    else:
        profile_pipeline.correct(args.output_dir, cfg or {})
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
