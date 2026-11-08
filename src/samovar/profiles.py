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
    """Read canonical CSV/TSV, CAMI or native profiler output at one rank."""
    path = Path(path)
    sidecar = Path(str(path) + ".profile.json")
    saved = json.loads(sidecar.read_text()) if sidecar.is_file() else {}
    meta = metadata({**saved, **dict(settings or {})})
    fmt = meta.get("format", "auto")
    with path.open() as handle:
        lines = handle.readlines()
    if fmt == "auto" and any(x.startswith("#consensus_taxonomy\tNCBI_tax_id\t") for x in lines):
        fmt = "motus"
    if fmt == "auto" and any('taxonomy_id' in x and 'rel_taxon_abundance' in x for x in lines[:1]):
        fmt = "phanta"
    if fmt == "auto" and any("taxonomy_id\ttaxonomy_lvl\t" in x and "new_est_reads" in x for x in lines):
        fmt = "bracken"
    if fmt == "auto" and any(x.startswith("Sample_file\t") and "Taxonomic_abundance" in x for x in lines):
        fmt = "sylph"
    if fmt == "auto":
        fmt = "cami" if any(x.startswith("@@TAXID") for x in lines) else (
            "metaphlan" if any("clade_name" in x and "relative_abundance" in x for x in lines) else "abundance")
    if fmt == "abundance":
        table = pd.read_csv(path, sep="\t" if path.suffix == ".tsv" else ",", dtype={"taxid": str})
    elif fmt == "phanta":
        if meta['taxonomy'] != 'ncbi' or meta['rank'] != 'species':
            raise ValueError('Phanta exposes species profiles with explicit unresolved native taxids')
        from samovar.taxdump import names_dmp, load_scientific_names
        names = names_dmp()
        canonical = load_scientific_names(names) if names is not None else {}
        required = ['name', 'taxonomy_id', 'new_est_reads', 'rel_taxon_abundance']
        native = pd.read_csv(path, sep='\t', dtype={'taxonomy_id': str}) if lines else pd.DataFrame(columns=required)
        if not set(required) <= set(native):
            raise ValueError('Malformed Phanta scaled Bracken profile')
        if native['taxonomy_id'].duplicated().any():
            raise ValueError('Duplicate Phanta species taxid')
        abundance = pd.to_numeric(native['rel_taxon_abundance'], errors='raise').astype(float)
        counts = pd.to_numeric(native['new_est_reads'], errors='raise').astype(float)
        if not np.isfinite(counts).all() or (counts < 0).any():
            raise ValueError('Invalid native Phanta read counts')
        ids, original = [], {}
        for row in native.itertuples(index=False):
            tax, name = str(row.taxonomy_id), str(row.name)
            if not tax.isdigit() or int(tax) < 1:
                raise ValueError('Invalid native Phanta species taxid')
            mapped = tax if canonical.get(int(tax)) == name else 'unmapped:phanta:' + tax
            ids.append(mapped)
            original[mapped] = dict(taxid=tax, name=name)
        sample = str(meta.get('sample', 'sample'))
        table = pd.DataFrame({'taxid': ids, 'N_' + sample: abundance})
        meta.update(unit='fraction', basis='organisms', normalization='classified', native_taxa=original)
    elif fmt in {"motus", "motus3"}:
        if meta["taxonomy"] != "ncbi":
            raise ValueError("mOTUs3 profiles use NCBI taxonomy with explicit unmapped mOTUs")
        import re
        header = next((line.rstrip("\n").split("\t") for line in lines
                       if line.startswith("#consensus_taxonomy\tNCBI_tax_id\t")), None)
        if header is None or len(header) != 3:
            raise ValueError("mOTUs3 requires a single-sample mOTU profile with -p NCBI taxids")
        if any(re.search(r"(?:^|\s)-c(?:\s|$)", line) for line in lines if line.startswith("#")):
            raise ValueError("mOTUs3 count output (-c) is not organism relative abundance")
        sample = str(meta.get("sample", header[2]))
        rows, seen = [], set()
        for line in lines:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 3:
                raise ValueError("Malformed mOTUs3 profile row")
            label, ncbi, value = fields
            match = re.search(r"\[((?:ref|meta|ext)_mOTU_[^]\s]+)\]$", label)
            native_id = match[1] if match else label
            if native_id in seen:
                raise ValueError(f"Duplicate mOTUs3 identifier: {native_id}")
            seen.add(native_id)
            if label == "unassigned":
                taxid = "0"
            elif not match:
                raise ValueError("mOTUs3 profile must retain mOTU IDs (-k mOTU)")
            elif native_id.startswith("ref_") and ncbi.isdigit() and int(ncbi) > 0:
                taxid = str(int(ncbi))
            else:
                # Meta/ext-mOTUs may name only an ancestral NCBI taxon, not a species.
                taxid = "unmapped:" + native_id
            rows.append((taxid, float(value)))
        col = "N_" + sample
        table = pd.DataFrame(rows, columns=["taxid", col])
        meta.update(kind="abundance", unit="fraction", basis="organisms", taxonomy="ncbi",
                    rank="species", normalization="all", native_header=[x.strip() for x in lines if x.startswith("#")])
        # Several reference mOTUs can resolve to the same NCBI species.
        values = pd.to_numeric(table[col], errors="raise")
        if not np.isfinite(values).all() or (values < 0).any():
            raise ValueError("Profile values must be finite and nonnegative")
        table = table.groupby("taxid", as_index=False, sort=False)[col].sum()
        meta["unassigned_by_sample"] = {sample: float(table.loc[table.taxid.eq("0"), col].sum())}
        table = table.loc[table[col] > 0].reset_index(drop=True)
    elif fmt == "bracken":
        raw = pd.read_csv(path, sep="\t", dtype={"taxonomy_id": str})
        required = {"taxonomy_id", "taxonomy_lvl", "new_est_reads"}
        if not required <= set(raw) or raw[list(required)].isna().any().any():
            raise ValueError("Expected Bracken abundance TSV with taxonomy_id, taxonomy_lvl, new_est_reads")
        if not raw.taxonomy_id.map(lambda value: value.isdigit() and int(value) > 0).all():
            raise ValueError("Bracken requires positive NCBI taxids")
        if not raw.taxonomy_lvl.eq("S").all():
            raise ValueError("Bracken profiler currently requires species (-l S)")
        table = raw[["taxonomy_id", "new_est_reads"]].rename(columns={
            "taxonomy_id": "taxid", "new_est_reads": "N_" + str(meta.get("sample", path.stem))})
        meta.update(kind="abundance", unit="counts", basis="reads", taxonomy="ncbi", rank="species", normalization="classified")
    elif fmt == "sylph":
        if meta["taxonomy"] != "ncbi" or meta["kind"] != "abundance":
            raise ValueError("Sylph genome mapping requires an NCBI abundance profile")
        mapping_path = meta.get("genome_taxonomy")
        if not mapping_path:
            raise ValueError("Sylph needs an explicit genome_taxonomy TSV (Genome_file, taxid); no taxids are inferred")
        mapping = pd.read_csv(mapping_path, sep="\t", dtype=str)
        if not {"Genome_file", "taxid"} <= set(mapping) or mapping[["Genome_file", "taxid"]].isna().any().any() or mapping.Genome_file.duplicated().any():
            raise ValueError("Invalid Sylph genome taxonomy mapping; require unique Genome_file and taxid")
        if not mapping.taxid.map(lambda value: value.isdigit() and int(value) > 0).all():
            raise ValueError("Sylph mapping requires positive NCBI taxids")
        sample = str(meta.get("sample", path.stem))
        if not lines:
            table = pd.DataFrame(columns=["taxid", "N_" + sample])
        else:
            raw = pd.read_csv(path, sep="\t", dtype={"Genome_file": str})
            if not {"Genome_file", "Sample_file", "Taxonomic_abundance"} <= set(raw):
                raise ValueError("Expected sylph profile output, not query output")
            if raw.Sample_file.isna().any() or raw.Sample_file.nunique() > 1:
                raise ValueError("Profiler execution must produce one Sylph sample")
            raw["taxid"] = raw.Genome_file.map(mapping.set_index("Genome_file").taxid)
            if raw.taxid.isna().any():
                missing = raw.loc[raw.taxid.isna(), "Genome_file"].tolist()
                raise ValueError(f"Sylph genomes lack explicit taxonomy mapping: {missing}")
            raw["Taxonomic_abundance"] = pd.to_numeric(raw.Taxonomic_abundance, errors="raise")
            if not np.isfinite(raw.Taxonomic_abundance).all() or (raw.Taxonomic_abundance < 0).any():
                raise ValueError("Sylph abundances must be finite and nonnegative")
            # Native percentages are printed to four decimals; remove only a
            # rounding-sized surplus before the common mass validation.
            total = raw.Taxonomic_abundance.sum()
            if 100 < total <= 100 + len(raw) * 0.00005 + 1e-9:
                raw["Taxonomic_abundance"] *= 100 / total
            table = raw.groupby("taxid", as_index=False).Taxonomic_abundance.sum().rename(columns={"Taxonomic_abundance": "N_" + sample})
        # A genome's NCBI ID can be below species. Existing quantitative stages
        # project taxid-rank inputs to species through the common taxonomy code.
        meta.update(unit="percent", basis="organisms", normalization="classified", rank="taxid")
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


def _fasta_length(path):
    from samovar.seqio import fasta_handle
    with fasta_handle(path) as handle:
        return sum(len(line.strip()) for line in handle if not line.startswith(">"))


def _reference_key(taxid):
    return str(taxid).split(".")[0]


def _explicit_reference_paths(taxids, config):
    """``profile_references`` maps one taxid to one genome. Length is not stored."""
    chosen = (config or {}).get("profile_references")
    if not chosen:
        return None
    if not isinstance(chosen, dict):
        raise ValueError("profile_references must map a taxid to one genome FASTA")
    wanted = {_reference_key(tax) for tax in taxids}
    paths = {}
    for tax in sorted(wanted):
        raw = chosen.get(tax)
        if isinstance(raw, (list, tuple, dict)):
            raise ValueError(
                f"Ambiguous profile reference for taxid {tax}: name one genome. "
                "Its length is read from the file."
            )
        if not raw:
            continue
        path = Path(str(raw)).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Profile reference genome for taxid {tax} is not a file: {path}")
        paths[tax] = str(path)
    missing = wanted - set(paths)
    if missing:
        raise ValueError(
            "profile_references has no genome for taxa: "
            + ", ".join(sorted(missing))
            + ". Name one FASTA per taxid; its length is read from the file."
        )
    return paths


def _local_reference_candidates(genome_dir, taxid):
    root = Path(genome_dir)
    if not root.is_dir():
        return []
    suffixes = (".fna.gz", ".fa.gz", ".fasta.gz", ".fsa.gz", ".fna", ".fa", ".fasta", ".fsa")
    hits = []
    for path in sorted(root.iterdir()):
        if not path.is_file() or path.name.startswith("."):
            continue
        name = path.name
        stem = next((name[: -len(ext)] for ext in suffixes if name.endswith(ext)), "")
        if stem == taxid or stem.startswith(taxid + "_"):
            hits.append(path)
    return hits


def _selected_reference_paths(taxids, genome_dir, saved):
    """One local FASTA, or the pinned file from an earlier explicit choice."""
    paths = {}
    for tax in sorted({_reference_key(tax) for tax in taxids}):
        hits = _local_reference_candidates(genome_dir, tax)
        if len(hits) > 1:
            listed = ", ".join(str(path) for path in hits)
            raise ValueError(
                f"Ambiguous reference genomes for taxid {tax}: {listed}. "
                "Set profile_references to one FASTA; its length is read from the file."
            )
        if len(hits) == 1:
            paths[tax] = str(hits[0].resolve())
            continue
        record = saved.get(tax)
        if record and record.get("path"):
            paths[tax] = str(Path(record["path"]).resolve())
    return paths


def reference_lengths(taxids, genome_dir, config=None, *, selected_paths=None):
    """Lengths are the base counts of the selected reference FASTA files.

    The public choice is ``profile_references`` (taxid → one genome) or the
    genomes already selected for this run. A directory with several FASTA
    files for one taxid is an error, not a silent cache pick.
    """
    genome_dir = Path(genome_dir)
    genome_dir.mkdir(parents=True, exist_ok=True)
    manifest = genome_dir / ".profile_references.json"
    saved = json.loads(manifest.read_text()) if manifest.is_file() else {}
    wanted = {_reference_key(tax) for tax in taxids}
    if selected_paths is not None:
        paths = {
            _reference_key(tax): str(Path(selected_paths[key]).expanduser().resolve())
            for tax in wanted
            for key in (tax, tax.split(".")[0])
            if key in selected_paths
        }
    else:
        paths = _explicit_reference_paths(wanted, config)
        if paths is None:
            paths = _selected_reference_paths(wanted, genome_dir, saved)
    sizes = {}
    for tax, raw in paths.items():
        path = Path(raw)
        if not path.is_file():
            raise FileNotFoundError(f"Profile reference genome for taxid {tax} is not a file: {path}")
        record = saved.get(tax) or {}
        same = record.get("path") and Path(record["path"]).resolve() == path.resolve()
        stat = path.stat()
        if same and (stat.st_size != record.get("size") or stat.st_mtime_ns != record.get("mtime_ns")):
            raise ValueError(
                f"Selected profile reference changed for taxid {tax}: {path}; rebuild dependent profiles"
            )
        length = _fasta_length(path)
        if same and record.get("length") not in (None, length):
            raise ValueError(
                f"Selected profile reference changed for taxid {tax}: {path}; rebuild dependent profiles"
            )
        sizes[tax] = length
        saved[tax] = dict(path=str(path.resolve()), length=length, size=stat.st_size, mtime_ns=stat.st_mtime_ns)
    missing = wanted - set(sizes)
    if missing:
        raise ValueError(
            "No explicit reference genome for profile taxa: "
            + ", ".join(sorted(missing))
            + ". Set profile_references to one FASTA per taxid; its length is read from the file."
        )
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
    def descriptor(value, label):
        if not isinstance(value, dict) or not value.get("path"):
            raise ValueError(f"{label} requires a path and explicit metadata")
        result = dict(value)
        result["path"] = str(Path(result["path"]).expanduser().resolve())
        for key in ("kind", "unit", "basis", "taxonomy", "rank", "normalization"):
            if key not in result:
                raise ValueError(f"{label} requires explicit {key}")
        metadata(result)
        if not Path(result["path"]).is_file():
            raise ValueError(f"{label}: file does not exist: {result['path']}")
        load_profile(result["path"], result)
        return result
    inputs = {str(name): descriptor(value, f"profile_inputs {name}") for name, value in inputs.items()}
    truth = descriptor(config["profile_truth"], "profile_truth") if config.get("profile_truth") else None
    corrector = config.get("profile_corrector", "profile_ratio")
    if corrector not in {"profile_ratio", "none"}:
        from samovar.abundance_correctors import require_known_export
        corrector = require_known_export(corrector)
    references = config.get("profile_references") or {}
    if references and not isinstance(references, dict):
        raise ValueError("profile_references must map a taxid to one genome FASTA")
    resolved = {}
    for tax, raw in references.items():
        key = _reference_key(tax)
        if isinstance(raw, (list, tuple, dict)) or not raw:
            raise ValueError(
                f"Ambiguous profile reference for taxid {key}: name one genome. "
                "Its length is read from the file."
            )
        path = Path(str(raw)).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"Profile reference genome for taxid {key} is not a file: {path}")
        if key in resolved:
            raise ValueError(f"Ambiguous profile reference for taxid {key}: name one genome.")
        resolved[key] = str(path)
    return inputs, truth, corrector, resolved


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


COMPOSITION_TOP_N = 10
COMPOSITION_COLORS = (
    "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
    "#a6761d", "#1f78b4", "#666666", "#e41a1c", "#377eb8", "#4daf4a",
    "#984ea3", "#ff7f00", "#a65628", "#f781bf", "#999999", "#66c2a5",
    "#fc8d62", "#8da0cb",
)


def composition_color(taxid: str) -> str:
    """Same taxon keeps the same colour across stages and annotators."""
    token = str(taxid)
    if token in {"0", "0 (unclassified)"}:
        return "#bdbdbd"
    if token == "Other":
        return "#eeeeee"
    if token.startswith("unmapped:") or token.endswith("(unresolved)"):
        return "#636363"
    import hashlib
    slot = int(hashlib.md5(token.encode()).hexdigest()[:8], 16) % len(COMPOSITION_COLORS)
    return COMPOSITION_COLORS[slot]


def composition_label(taxid: str) -> str:
    token = str(taxid)
    if token == "0":
        return "0 (unclassified)"
    if token.startswith("unmapped:"):
        return f"{token} (unresolved)"
    return token


def composition_fractions(table, info):
    """Column shares on the table's own scale. Unclassified stays taxid 0.

    A classified profile that omitted taxid 0 keeps ``unassigned_by_sample`` as
    that slice and scales the other taxa only when they already sum to one.
    Returns None for presence profiles, which have no quantitative composition.
    """
    if info.get("kind") != "abundance":
        return None
    work, info = validate(table, info)
    cols = n_sample_columns(work)
    values = work[cols]
    if info["unit"] == "percent":
        work[cols] = values / 100
    elif info["unit"] == "counts":
        work[cols] = values.div(values.sum().replace(0, np.nan)).fillna(0)
    samples = [column[2:] for column in cols]
    unassigned = {
        sample: float((info.get("unassigned_by_sample") or {}).get(sample) or 0)
        for sample in samples
    }
    if not work.taxid.eq("0").any() and any(unassigned.values()):
        for column, sample in zip(cols, samples):
            total = float(work[column].sum())
            room = 1 - unassigned[sample]
            if total > room + 1e-6 and total > 0 and room >= 0:
                work[column] = work[column] * (room / total)
        extra = {"taxid": "0"}
        extra.update({column: unassigned[sample] for column, sample in zip(cols, samples)})
        work = pd.concat([work, pd.DataFrame([extra])], ignore_index=True)
    return work


def composition_display(fractions_table, top_n=COMPOSITION_TOP_N):
    """Top taxa plus Other. Unclassified and unresolved stay out of Other."""
    cols = n_sample_columns(fractions_table)
    labels = fractions_table.taxid.astype(str)
    special = labels.eq("0") | labels.str.startswith("unmapped:")
    ranked = fractions_table.loc[~special].copy()
    ranked["_mean"] = ranked[cols].mean(axis=1) if cols else 0
    ranked = ranked.sort_values(["_mean", "taxid"], ascending=[False, True])
    keep = ranked.head(int(top_n)).drop(columns="_mean")
    rest = ranked.iloc[int(top_n):]
    shown = pd.concat([keep, fractions_table.loc[special]], ignore_index=True)
    if not rest.empty:
        other = {"taxid": "Other"}
        for column in cols:
            other[column] = float(rest[column].sum())
        shown = pd.concat([shown, pd.DataFrame([other])], ignore_index=True)
    return shown


def alpha_metrics(table, info):
    """Observed richness and Shannon for one profile, one sample at a time.

    Shannon is the existing ``shannon_entropy`` (natural log, nats) on positive
    biological taxa. Taxid 0 and ``unmapped:`` stay out of both metrics.
    Presence profiles have no abundance and return None.
    """
    full = composition_fractions(table, info)
    if full is None:
        return None
    from samovar.table_scorers import shannon_entropy

    labels = full.taxid.astype(str)
    biological = ~(labels.eq("0") | labels.str.startswith("unmapped:"))
    rows = []
    for column in n_sample_columns(full):
        bio = full.loc[biological, column].to_numpy(dtype=float)
        positive = bio[np.isfinite(bio) & (bio > 0)]
        rows.append({
            "sample": column[2:],
            "richness": int(positive.size),
            "shannon": float(shannon_entropy(positive)),
            "unclassified": float(full.loc[labels.eq("0"), column].sum()) if labels.eq("0").any() else None,
            "unresolved": float(full.loc[labels.str.startswith("unmapped:"), column].sum()) if labels.str.startswith("unmapped:").any() else 0.0,
        })
    return rows


BETA_METRICS = ("bray_curtis",)
BETA_PREPROCESSING = (
    "full column fractions; every taxon, including taxid 0 and unmapped; "
    "zero-mass samples omitted rather than given distance 0; no top-N; no rarefaction"
)


def classical_pcoa(distance):
    """Classical PCoA. Negative eigenvalues stay negative and are not axes."""
    matrix = np.asarray(distance, dtype=float)
    size = matrix.shape[0]
    centered = -0.5 * matrix ** 2
    scale = np.full((size, size), 1.0 / size)
    gram = centered - scale @ centered - centered @ scale + scale @ centered @ scale
    values, vectors = np.linalg.eigh(gram)
    order = np.argsort(values)[::-1]
    values = values[order]
    vectors = vectors[:, order]
    scale = max(1e-10, 1e-8 * float(np.max(np.abs(values))) if values.size else 1e-10)
    positive = values > scale
    negative = values < -scale
    coordinates = {}
    explained = []
    positive_sum = float(values[positive].sum()) if positive.any() else 0.0
    for axis, index in enumerate(np.where(positive)[0], start=1):
        coordinates[f"PC{axis}"] = (vectors[:, index] * np.sqrt(values[index])).tolist()
        explained.append(float(values[index] / positive_sum) if positive_sum else 0.0)
    return {
        "eigenvalues": [float(value) for value in values],
        "negative_eigenvalues": [float(value) for value in values[negative]],
        "explained_variance": explained,
        "coordinates": coordinates,
    }


def beta_diversity(table, info, metric="bray_curtis"):
    """Distance matrix and PCoA for one profile. Samples are not pooled across profiles."""
    if metric not in BETA_METRICS:
        raise ValueError(f"Unknown beta metric {metric!r}. Known metrics: {', '.join(BETA_METRICS)}")
    result = {
        "metric": metric,
        "preprocessing": BETA_PREPROCESSING,
        "rank": info.get("rank", ""),
        "taxonomy": info.get("taxonomy", ""),
        "unit": info.get("unit", ""),
        "basis": info.get("basis", ""),
        "normalization": info.get("normalization", ""),
        "sample_ids": [],
        "omitted_sample_ids": [],
        "matrix": None,
        "eigenvalues": [],
        "negative_eigenvalues": [],
        "explained_variance": [],
        "coordinates": {},
        "status": None,
    }
    full = composition_fractions(table, info)
    if full is None:
        result["status"] = "presence"
        return result
    included, vectors = [], []
    for column in n_sample_columns(full):
        sample = column[2:]
        vector = full[column].to_numpy(dtype=float)
        if float(np.nansum(vector)) <= 0:
            result["omitted_sample_ids"].append(sample)
            continue
        included.append(sample)
        vectors.append(vector)
    result["sample_ids"] = included
    if len(included) < 2:
        result["status"] = "fewer than two samples with positive mass"
        return result
    matrix = _distance_matrix(vectors, metric)
    result["matrix"] = matrix.tolist()
    result.update(classical_pcoa(matrix))
    if float(matrix.max()) <= 1e-12 or not result["explained_variance"]:
        result["status"] = "degenerate distances"
    elif len(result["explained_variance"]) < 2:
        result["status"] = "one positive eigenvalue"
    return result


GROUPED_PREPROCESSING = (
    BETA_PREPROCESSING
    + "; one matrix for sample × profile points that share taxonomy, rank, basis, and normalization; "
    "common samples are the intersection of positive-mass samples; "
    "a sample missing from one profile is excluded, not filled with zeros; "
    "a taxon missing from one profile is zero after alignment"
)


def _positive_fraction_series(table, info):
    full = composition_fractions(table, info)
    if full is None:
        return None
    labels = full.taxid.astype(str)
    series = {}
    for column in n_sample_columns(full):
        values = pd.Series(full[column].to_numpy(dtype=float), index=labels)
        if float(values.sum()) <= 0:
            continue
        series[column[2:]] = values
    return series


def _distance_matrix(vectors, metric):
    if metric != "bray_curtis":
        raise ValueError(f"Unknown beta metric {metric!r}. Known metrics: {', '.join(BETA_METRICS)}")
    from samovar.sample_scorers import bray_curtis_distance

    size = len(vectors)
    matrix = np.zeros((size, size))
    for left in range(size):
        for right in range(left + 1, size):
            distance = float(bray_curtis_distance(vectors[left], vectors[right]))
            matrix[left, right] = matrix[right, left] = distance
    return matrix


def grouped_beta(tables, metas, metric="bray_curtis"):
    """One distance matrix and PCoA per compatible set of profiles.

    The group is the profile identity. Stages are not mixed here; pass one stage.
    """
    if metric not in BETA_METRICS:
        raise ValueError(f"Unknown beta metric {metric!r}. Known metrics: {', '.join(BETA_METRICS)}")
    classes = {}
    excluded = []
    for name, table in tables.items():
        info = dict(metas.get(name) or {})
        vectors = _positive_fraction_series(table, info)
        if vectors is None:
            excluded.append({"profile": str(name), "reason": "presence"})
            continue
        key = (
            str(info.get("taxonomy") or ""),
            str(info.get("rank") or ""),
            str(info.get("basis") or ""),
            str(info.get("normalization") or ""),
        )
        classes.setdefault(key, []).append((str(name), info, vectors))
    spaces = []
    for key, members in classes.items():
        common = set.intersection(*(set(vectors) for _, _, vectors in members))
        omitted = [
            {"profile": name, "sample": sample, "reason": "not in every compatible profile"}
            for name, _, vectors in members
            for sample in sorted(set(vectors) - common)
        ]
        points = []
        for name, info, vectors in sorted(members, key=lambda item: item[0]):
            for sample in sorted(common):
                points.append({
                    "sample": sample,
                    "profile": name,
                    "tool": str(info.get("tool") or ""),
                    "run_name": str(info.get("run_name") or ""),
                    "database": str(info.get("database") or ""),
                    "vector": vectors[sample],
                })
        space = {
            "metric": metric,
            "preprocessing": GROUPED_PREPROCESSING,
            "taxonomy": key[0],
            "rank": key[1],
            "basis": key[2],
            "normalization": key[3],
            "profiles": [name for name, _, _ in members],
            "excluded_profiles": list(excluded),
            "omitted_samples": omitted,
            "points": [],
            "matrix": None,
            "eigenvalues": [],
            "negative_eigenvalues": [],
            "explained_variance": [],
            "coordinates": {},
            "status": None,
        }
        if len(points) < 2:
            space["points"] = [{key: point[key] for key in ("sample", "profile", "tool", "run_name", "database")} for point in points]
            space["status"] = "fewer than two points"
            spaces.append(space)
            continue
        taxa = sorted(set().union(*(set(point["vector"].index) for point in points)))
        aligned = [point["vector"].reindex(taxa, fill_value=0.0).to_numpy(dtype=float) for point in points]
        matrix = _distance_matrix(aligned, metric)
        space["matrix"] = matrix.tolist()
        space.update(classical_pcoa(matrix))
        if float(matrix.max()) <= 1e-12 or not space["explained_variance"]:
            space["status"] = "degenerate distances"
        elif len(space["explained_variance"]) < 2:
            space["status"] = "one positive eigenvalue"
        exported = []
        for index, point in enumerate(points):
            row = {key: point[key] for key in ("sample", "profile", "tool", "run_name", "database")}
            for axis, values in space["coordinates"].items():
                row[axis] = values[index]
            exported.append(row)
        space["points"] = exported
        spaces.append(space)
    return {"spaces": spaces, "excluded_profiles": excluded}


def compatible(pred, truth):
    for field in ("taxonomy", "rank"):
        if pred[field] != truth[field]:
            raise ValueError(f"Incompatible profile {field}: {pred[field]} vs {truth[field]}")
    if truth["kind"] != "presence":
        for field in ("basis", "normalization"):
            if pred[field] != truth[field]:
                raise ValueError(f"Incompatible profile {field}: {pred[field]} vs {truth[field]}")


def report_profiles(annotation, dest, stage, config=None, *, annotators="all"):
    """Publish the report's typed profile contract and return per-sample inventory.

    No rank/basis conversion, pooling, GT lookup or top-N filtering belongs at
    this boundary. Binner identities retain the full source column so different
    database runs cannot overwrite one another through display-name parsing.
    """
    from samovar.abundance import annotation_to_abundance_tables, taxid_annotation_columns, annotator_name

    if stage not in {"initial", "regenerated", "reprofiled"}:
        raise ValueError(f"Invalid report profile stage: {stage}")
    tables, metas = {}, {}
    frame = annotation.DataFrame
    cached_names = {annotator_name(column) for column in taxid_annotation_columns(frame)}
    for column in taxid_annotation_columns(frame):
        columns = [column] + (["sample"] if "sample" in frame else [])
        derived = annotation_to_abundance_tables(frame[columns])
        if not derived:
            continue
        name = str(column)
        tables[name] = next(iter(derived.values()))
        metas[name] = metadata(dict(rank="taxid", normalization="all",
                                    representation="per-read", source_column=name))
    for name, table in (annotation.abundance_tables or {}).items():
        if name not in annotation.profile_metadata and name in cached_names:
            continue  # abundance cache of the read columns above
        if name in tables:
            raise ValueError(f"Report profile identity collision: {name}")
        tables[name] = table
        metas[name] = metadata(annotation.profile_metadata.get(name, {"rank": "taxid"}))
        metas[name]["representation"] = "profile"
    runs = (config or {}).get("run_config") or []
    rows = {}
    selected_tables, selected_metas = {}, {}
    for name, table in tables.items():
        info = metas[name]
        source = info.get("source_identity", name)
        matches = [run for run in runs if source == run.get("run_name") or
                   source.startswith("taxID_" + str(run.get("run_name")) + "_")]
        run = matches[0] if len(matches) == 1 else {}
        info.update(stage=stage, identity=name, run_name=info.get("run_name", run.get("run_name", source)),
                    tool=info.get("tool", run.get("type")),
                    database=info.get("database", run.get("db_path")))
        from samovar.multiqc_report import profile_selected
        if not profile_selected(name, info, annotators):
            continue
        table, info = validate(table, info)
        # Metadata unassigned values are fractions, independently of table units.
        unassigned = dict(info.get("unassigned_by_sample") or {})
        samples = [c[2:] for c in n_sample_columns(table)]
        if set(unassigned) - set(samples):
            raise ValueError(f"Unassigned metadata has unknown samples for {name}")
        relative = fractions(table, {**info, "normalization": "all"}) if info["kind"] == "abundance" else None
        for col, sample in zip(n_sample_columns(table), samples):
            mass = float(table[col].sum())
            value = unassigned.get(sample)
            if relative is not None and table.taxid.eq("0").any():
                observed = float(relative.loc[relative.taxid == "0", col].sum())
                if value is not None and not np.isclose(value, observed):
                    raise ValueError(f"Conflicting unassigned mass for {name}/{sample}")
                value = observed
            if value is not None:
                if not np.isfinite(value) or not 0 <= value <= 1:
                    raise ValueError(f"Invalid unassigned fraction for {name}/{sample}")
                unassigned[sample] = float(value)
            unresolved = float(relative.loc[table.taxid.str.startswith("unmapped:"), col].sum()) if relative is not None else None
            status = "all-unclassified" if value == 1 else ("empty" if mass == 0 else "available")
            rows[str(len(rows))] = dict(profile=name, sample=sample, stage=stage,
                tool=info["tool"], run_name=info["run_name"], database=info["database"],
                method=info.get("method"), representation=info["representation"],
                kind=info["kind"], unit=info["unit"], basis=info["basis"], taxonomy=info["taxonomy"],
                rank=info["rank"], normalization=info["normalization"], mass=mass,
                unassigned_fraction=value, unresolved_fraction=unresolved, status=status)
        info["unassigned_by_sample"] = unassigned
        selected_tables[name], selected_metas[name] = table, info
    write_profiles(dest, selected_tables, selected_metas)
    return rows



def run_motus(r1, r2, database, raw_output, sample, threads=1, executable="motus", flags=""):
    """Run mOTUs3 with its configured companions and a private native workspace."""
    import argparse
    import os
    import shlex
    import subprocess
    import tempfile
    from samovar.paths import command_argv, resolve_executable
    from samovar.seqio import fastq_has_reads, is_gzip_file

    options = argparse.ArgumentParser(add_help=False)
    options.add_argument("-g", type=int, default=3)
    options.add_argument("-l", type=int, default=75)
    options.add_argument("-v", type=int, default=3)
    opts = options.parse_args(shlex.split(flags))
    if opts.g < 1 or opts.l < 1 or threads < 1:
        raise ValueError("mOTUs3 requires positive marker, alignment-length and thread settings")
    db, out = Path(database).resolve(), Path(raw_output).resolve()
    versions = list(db.glob("*_versions"))
    if not versions and (db / "db_mOTU").is_dir():
        db = db / "db_mOTU"
        versions = list(db.glob("*_versions"))
    if len(versions) != 1:
        raise FileNotFoundError(f"mOTUs3 database requires one *_versions file in {db}")
    prefix = versions[0].name[:-len("_versions")]
    info = dict(line.split("\t", 1) for line in versions[0].read_text().splitlines() if "\t" in line)
    if not info.get("motus", "").startswith("3."):
        raise ValueError("mOTUs3 requires a version 3 database")
    if not (db / (prefix + "_DB_NR.fasta.bwt")).is_file():
        raise FileNotFoundError("mOTUs3 database is missing its BWA marker index")
    command = command_argv(resolve_executable(executable, tool_key="motus"))
    companions = {name: command_argv(resolve_executable(name, tool_key=name))
                  for name in ("bwa", "samtools")}
    for argv in (command, *companions.values()):
        if not argv or not Path(argv[0]).is_file():
            raise FileNotFoundError(f"mOTUs3 executable is not configured: {argv}")
    out.parent.mkdir(parents=True, exist_ok=True)
    scratch = out.parent / ".tmp" / "motus"
    scratch.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=scratch, prefix="profile-") as tmp:
        work = Path(tmp)
        # Native mOTUs3 hardcodes bwa/samtools names and splits paths internally.
        # Expose only config-resolved commands, and pass safe relative input paths.
        bindir = work / "bin"; bindir.mkdir()
        for name, argv in companions.items():
            launcher = bindir / name
            launcher.write_text("#!/bin/sh\nexec " + shlex.join(argv) + ' "$@"\n')
            launcher.chmod(0o700)
        env = {**os.environ, "PATH": str(bindir) + os.pathsep + os.defpath}
        (work / prefix).symlink_to(db, target_is_directory=True)
        inputs = []
        for number, source in enumerate([r1, *([r2] if r2 else [])]):
            filename = f"reads{number}.fastq" + (".gz" if is_gzip_file(source) else "")
            (work / filename).symlink_to(Path(source).resolve())
            inputs.extend((["-f", filename] if number == 0 and r2 else ["-r" if number else "-s", filename]))
        target = work / "profile.tsv"
        if fastq_has_reads(r1):
            argv = [*command, "profile", *inputs, "-db", prefix, "-n", sample, "-p", "-k", "mOTU",
                    "-t", str(threads), "-g", str(opts.g), "-l", str(opts.l), "-v", str(opts.v), "-o", target.name]
            subprocess.run(argv, cwd=work, env=env, check=True)
            # An rc=0 without a valid typed result is still an error.
            load_profile(target, {"format": "motus", "sample": sample})
        else:
            target.write_text(f"#consensus_taxonomy\tNCBI_tax_id\t{sample}\n")
        target.replace(out)


def run_bracken(report, database, raw_output, read_length=100, threshold=10):
    """Run the native estimator via configured Python; empty eligible sets are valid."""
    import subprocess
    from samovar.paths import command_argv, python_path, resolve_executable
    distribution = Path(database) / f"database{read_length}mers.kmer_distrib"
    if not distribution.is_file() or not distribution.stat().st_size:
        raise FileNotFoundError(distribution)
    eligible = False
    for line in Path(report).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 6 or not fields[1].strip().isdigit():
            raise ValueError("Malformed Kraken2 report for Bracken")
        if fields[3] == "S" and int(fields[1]) > 0 and int(fields[1]) >= threshold:
            eligible = True
    if not eligible:
        Path(raw_output).write_text("taxonomy_id\ttaxonomy_lvl\tnew_est_reads\n")
        return
    script = resolve_executable("est_abundance.py", tool_key="est_abundance.py")
    subprocess.run([python_path(), *command_argv(script), "-i", str(report), "-o", str(raw_output),
                    "-k", str(distribution), "-l", "S", "-t", str(threshold)], check=True)


def main(argv=None):
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["ingest", "bracken", "motus", "collect", "materialize", "correct"])
    parser.add_argument("--source")
    parser.add_argument("--output")
    parser.add_argument("--sample")
    parser.add_argument("--identity")
    parser.add_argument("--database", default="")
    parser.add_argument("--genome-taxonomy")
    parser.add_argument("--kraken-report")
    parser.add_argument("--read-length", type=int, default=100)
    parser.add_argument("--threshold", type=int, default=10)
    parser.add_argument("--format", default="auto")
    parser.add_argument("--output-dir")
    parser.add_argument("--stage", default="initial")
    parser.add_argument("--config")
    parser.add_argument("--input-r1")
    parser.add_argument("--input-r2")
    parser.add_argument("--executable", default="motus")
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--native-flags", default="")
    args = parser.parse_args(argv)
    if args.action == "motus":
        run_motus(args.input_r1, args.input_r2, args.database, args.source, args.sample,
                  args.threads, args.executable, args.native_flags)
        return 0
    if args.action == "bracken":
        run_bracken(args.kraken_report, args.database, args.source, args.read_length, args.threshold)
        args.format = "bracken"
    if args.action in {"ingest", "bracken"}:
        if args.format == "abundance" and not Path(str(args.source) + ".profile.json").is_file():
            raise ValueError("Profiler contract requires an output .profile.json sidecar")
        settings = {"format": args.format, "sample": args.sample}
        if args.genome_taxonomy:
            settings["genome_taxonomy"] = args.genome_taxonomy
        table, meta = load_profile(args.source, settings)
        cols = n_sample_columns(table)
        if len(cols) != 1:
            raise ValueError("Profiler execution must produce one sample")
        table = table.rename(columns={cols[0]: "N_" + args.sample})
        meta.update(identity=args.identity, database=args.database, source=str(args.source))
        if args.action == "bracken":
            meta.update(read_length=args.read_length, threshold=args.threshold)
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
