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

Translate choices into the public prepare CLI: supplied reads use `--input_dir`;
generated data reuse their known output directory and generation config. Set
`--cores`, regeneration depth/sample count, selected annotators and any stage
endpoint from the agreed request. Preserve the selected read generator/type from
generation, including ONT/hybrid. Obtain exact flag spellings from current help.
Record input audit, scientific choices, argv and any warnings in `audit.md`.

A successful prepare produces configuration and `.log/samovar.sh` in the run;
it does not mean analysis has run. Inspect the generated config/script and
citations for the intended input, tools, database identities, resources and
stages. Do not suppress errors or edit generated configs to hide missing tools.
Continue with [exec](../samovar-exec/SKILL.md) only when execution was requested.
For an existing prepared run, compare its config/checkpoints first; preserve
completed results and request a branch/new output if the analysis changes.
