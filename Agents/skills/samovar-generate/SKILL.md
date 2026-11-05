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

Check sample labels against the installed workflow too: the current ISS workflow
expects `N_1` through `N_<n_samples>`. A table such as `N_0,N_1` can be accepted
at prepare time but fail during execution because generated and expected FASTQ
names differ. Warn and stop for a supported-input decision in that case; never
silently rename the user's samples. This is a known CLI limitation, not evidence
of successful generation.

Translate sample count, sequencing records, host fraction, seed and core budget
explicitly. ISS counts paired fragments; ONT counts single reads. Do not change
the simulator or drop host/taxa to make an unavailable mode work. Display and
record the exact argv and its interpretation in the run audit.

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
