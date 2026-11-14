---
name: samovar-prepare
description: Translate an explicit SAMOVAR analysis request into a prepared pipeline using verified reads, registered annotators and databases, and recorded scientific choices.
---

# Prepare a SAMOVAR run

Apply only to SAMOVAR workflow preparation. Read
[scientific integrity](references/rules/scientific-integrity-samovar.md) and
[code reuse](references/rules/code-usage.md).

Inspect `samovar prepare --help` and `samovar tools --status`. Establish whether
inputs are supplied reads or a previously prepared generated community. Require
the actual input/output directories and check read layout, sample identities,
annotators, database indexes and available resources. If real versus generated
inputs or catalog entries are ambiguous, ask before preparing; do not infer an
unmentioned generation step or choose an arbitrary database.

Named annotators are dynamic CLI options and may not appear in static help:
the form is `--ALIAS "TOOL INDEX"`, with a unique output alias and registered
tool/index names (for example `--kraken2-phage "kraken2 phage_test"`). Reuse
catalog identities; the alias does not create or install a new tool.
Combining annotators into an ensemble does not itself request the named
`--reprofiler ensemble` method. Preserve the installed default reprofiler unless
the user/profile explicitly selects a different implementation.

Check input ownership before execution. In SAMOVAR 0.13.2+, `sort_reads` writes
an owned sorted file atomically over the run's link, preserving external FASTQ
targets (including hardlinks and read-only sources). Older installations can
rewrite inputs through symlinks: warn and arrange an owned working copy there
when originals must remain unchanged. For large data, explain disk/copy costs
before copying blindly. Sorting still needs writable space inside the run, and
the first sorting stage can change hashes/mtimes of its owned output files.

Reuse registered tools and databases. If something required is missing, explain
what it is and use the existing installation/import dialogue only with the user's
authorization. Never replace an unavailable assembler or simulator with a different
analysis. Where no assembly is requested/available, inspect supported table
regeneration modes: `concat_weighted` pools abundance votes while `direct` keeps
per-annotator regeneration. If both are available and no approved policy or
explicit choice resolves them, ask; do not invent a scientific default. If only
one appropriate mode is available, explain it and confirm the analysis choice.

Ground truth has three equal inputs. Genome-based per-read taxids in FASTQ or FASTA headers are the default and preferred source. `--initial-ground-truth-table` is the explicit seq→taxid table and, with FASTQ, also a sample profile. `profile_truth` in the input YAML is the sample-level abundance or presence gold; it does not replace per-read labels. When several are set, profile scores stay separate and name `ground_truth`. `--repeat N` (CLI overrides YAML `repeat`) prepares each `out_XXXX` as its own run. Omit `--repeat` on an existing package to keep its N. A different N is refused. `--repeat N` on an existing FASTQ directory, including `1`, writes one with-replacement draw per child. Each draw has exactly as many single reads or paired fragments as the source. Paired mates share one index list and one new read id; ONT or other single-end files are drawn separately and are not mixed into the paired draw. Source FASTQ and source truth stay unchanged. A saved profile still cannot be repeated with N greater than 1; `--repeat 1` on a profile is one run, not a replicate series. A generated package is resumed, not nested. `--mode benchmark` prepares only the initial window, ending at `viz_initial` unless the request names an earlier endpoint. Do not add regeneration, training, or reprofiling to obtain an export. Each child receives that same core and memory budget and does not read another child's output. A directory path may contain spaces. A supplied-profile benchmark uses `profile_inputs` and `profile_truth` inside `--input_config`. Do not pass that per-read flag as a substitute for `profile_truth`.

Translate choices into the public prepare CLI: supplied reads use `--input_dir`;
generated data reuse their known output directory and generation config. Set
`--cores`, regeneration depth/sample count, selected annotators and any stage
endpoint from the agreed request. Preserve the selected read generator/type from
generation, including ONT/hybrid. Obtain exact flag spellings from current help.
Record input audit, scientific choices, argv and any warnings in `audit.md`.

Verify actual read counts rather than treating `N_reads` as a hard cap. With
`rescale_abundance=false`, tables retain their observed mass; pooled modes such
as `concat_weighted` can exceed that number. Current `concat_weighted` also
ignores `rescale_abundance=true`; setting it does not enforce a fixed depth.
If exact regenerated depth is required, establish whether the selected method
supports it and ask about an unsupported combination; do not silently substitute
another method. For methods that support rescaling, use the public
`--input_config` YAML interface and verify resulting tables/FASTQs. Do not patch
generated configuration or silently change abundance policy to pass a count check.

A successful prepare produces configuration and `.log/samovar.sh` in the run;
it does not mean analysis has run. Inspect the generated config/script and
citations for the intended input, tools, database identities, resources and
stages. Do not suppress errors or edit generated configs to hide missing tools.
Continue with [exec](../samovar-exec/SKILL.md) only when execution was requested.
For an existing prepared run, compare its config/checkpoints first; preserve
completed results and request a branch/new output if the analysis changes.
