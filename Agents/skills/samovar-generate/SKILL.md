---
name: samovar-generate
description: Prepare and execute synthetic metagenomic read generation with SAMOVAR from an explicit community or abundance table. Use for samovar generate requests.
---

# Generate synthetic reads with SAMOVAR

Apply only to SAMOVAR read/community generation. Read the bundled
[scientific integrity](references/rules/scientific-integrity-samovar.md) and
[code reuse](references/rules/code-usage.md) rules.

Inspect `samovar generate --help` and `samovar tools --status` for the installed
version. Require the community, simulator/read type and output directory. Resolve
taxon names to verified identifiers; ask if ambiguous. A request such as
“10+10+10” needs named groups and units before it can become abundances or read
counts. Never interpret missing quantities or an unspecified simulator silently.

Inspect supplied FASTA files or the abundance table. The current CLI accepts
`--genome_dir`, repeatable/comma-separated `--genomes` to lock FASTA paths/taxids,
and `--abundance` for a CSV with `taxid,N_<sample>…` columns. Use only options
shown by the installed help. For taxid requests, check local genomes/catalog
before any download; report unavailable taxa and obtain authorization for new
downloads. Bundled truncated genomes are test fixtures, not a scientific reference.

Check sample labels against the installed version too. From SAMOVAR 0.14.3,
ISS preserves literal `N_<sample>` labels, including zero-based and named samples,
and infers the count when `--n_samples` is omitted. An explicit count must match.
IDs start with a letter/digit and contain only letters/digits, `_` or `-`;
duplicates and unsafe IDs are errors. Older ISS workflows require `N_1` through
`N_<n_samples>` and can fail on other labels. Warn and stop on that older-version
limitation; never silently rename the user's samples or claim a failed run passed.

Translate sample count, sequencing records, host fraction, seed and core budget
explicitly. ISS counts paired fragments; ONT counts single reads. Do not change
the simulator or drop host/taxa to make an unavailable mode work. Display and
record the exact argv and its interpretation in the run audit.

`--repeat N` records a package and writes each `out_XXXX/.generate` through the existing generator. A supplied abundance table stays byte-for-byte that file; only the simulator seed changes. Omit `--repeat` for the single directory. Do not point `--repeat` at an existing single run or at a child `out_XXXX`. Executing the package script runs the child generators; it does not prepare the main pipeline.

`generate` prepares a generation workflow; it does not itself prove reads exist.
Execute the generated workflow through `samovar exec` only if authorized. A
generation-only directory can produce reads and then report that the main
pipeline has not been prepared. Inspect the actual status/artifacts and explain
this boundary; do not call the whole pipeline successful. For a requested full
analysis, continue with [prepare](../samovar-prepare/SKILL.md) before
[exec](../samovar-exec/SKILL.md).

Verify sample identities, actual nonempty FASTQ files, mate consistency, read
counts/types and requested truth outputs. Preserve provenance. For reuse, inspect
the existing generation config/checkpoints and compare the requested community
and parameters before resuming; do not overwrite a different run or use blanket
redo. Return paths to both the prepared workflow and the reads actually produced.
