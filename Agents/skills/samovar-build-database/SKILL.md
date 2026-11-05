---
name: samovar-build-database
description: Build or reuse a registered taxonomic annotator database with the SAMOVAR CLI, checking local genomes and taxonomy before downloads.
---

# Build a SAMOVAR annotator database

Apply only to SAMOVAR database selection/building. Read
[scientific integrity](references/rules/scientific-integrity-samovar.md) and
[code reuse](references/rules/code-usage.md).

Inspect `samovar build --help` and `samovar tools --status`. Require the annotator
type, verified genome/community source, database destination and, if registration
is requested, an unambiguous index name. Reject unknown types before launching
a builder. Current types include Kraken2, KrakenUniq, Kaiju and Centrifuge;
the installed help is authoritative.

Check for a compatible registered database first. Reuse it if it matches the
requested reference/taxonomy; an existing name alone does not prove compatibility.
Inspect the genome directories in the main install config and supplied build
YAML. Validate the requested taxon coverage against available FASTA/annotation
files. Reuse local valid genomes; use the public `samovar genome` interface for
authorized downloads if needed, never a private downloader or second cache.

Warn explicitly about missing taxa or failed downloads and stop for a decision
if this changes the requested database. Do not quietly build an incomplete index.
Keep a log of reused, downloaded and unavailable genomes. Do not synthesize taxids
or bypass taxonomy validation. Record builder errors, including OOM, as failures.

Use `samovar build --type … --config_path … --db_path …` with the supported
thread option. Inspect the build YAML schema in that installation's documentation
or prepared config; do not invent a schema. Add `--index` only when registering
the result is authorized. Do not overwrite another reference under that name.
Programs and taxonomy must resolve from SAMOVAR's install config.

After the build, check return status and type-specific index artifacts and
registration when requested. Record source genomes, taxonomy and omissions in
the run audit. On a second request, inspect and reuse a suitable completed index
instead of rebuilding/redownloading. Do not equate index files with validation
of every supported external builder; report which builder actually ran.
