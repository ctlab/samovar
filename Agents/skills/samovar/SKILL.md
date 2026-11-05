---
name: samovar
description: Plan and run metagenomic annotation or reprofiling with the SAMOVAR CLI, or inspect and resume an existing SAMOVAR run. Use for SAMOVAR analysis requests, not unrelated software development.
---

# SAMOVAR

Apply only when the user requests SAMOVAR analysis or help with an existing
SAMOVAR run. Read the bundled conditional rules before acting:

- [Scientific integrity](references/rules/scientific-integrity-samovar.md)
- [Reuse installed code and tools](references/rules/code-usage.md)

These files travel with the skill; the user does not install rules separately.

Establish the requested analysis, actual input files, output directory, required
annotators/databases, and resource limits from the request and existing run
configuration. Check files and installed tools; do not invent paths or taxonomic
identifiers. Ask about missing choices that change the scientific result. In
particular, synthetic generation needs an explicit community, simulator and
output location. Do not quietly select replacements for unavailable tools.

Use the installed `samovar help` and the relevant subcommand's `--help` to
construct commands for that installation. Inspect `samovar tools --status`
before installing or importing anything. If `samovar` is unavailable, report
that dependency and offer the repository's documented installation procedure.
Do not substitute inline Python or a private pipeline for a missing CLI stage.

Audit the supplied data and the request first, recording both in `audit.md` in
the agreed output directory. Route an authorized run through the installed
stage skills (they ship in the same bundle):

1. [Generate](../samovar-generate/SKILL.md) only for requested synthetic inputs;
   skip it when the user supplies reads.
2. [Build database](../samovar-build-database/SKILL.md) only if a suitable
   registered database is missing; otherwise record its reuse.
3. [Prepare](../samovar-prepare/SKILL.md) using the chosen data and tools.
4. [Exec](../samovar-exec/SKILL.md) through the requested endpoint, including
   checkup and checkpoint inspection before resume.

Follow these stage procedures in order; if the environment offers delegation
and the user authorizes it, they can be delegated with the same paths, choices
and audit. No IDE-specific subagent API is required. A prepared run's directory
must be known before execution. Do not delete checkpoints or request `--redo`
to hide a failure. If the request is ambiguous, ask a focused question and do
not create a run until the missing scientific choice is resolved.

Explain material parameter choices and record commands and decisions with the
run (an `audit.md` in its output directory is suitable). Execute only the work
the user authorized. Check return codes and requested artifacts, including
per-annotator outputs and regenerated data, before reporting completion. Report
partial, empty or failed stages explicitly; a zero exit code alone is insufficient.
Return artifact paths and the actual validation performed. Use the run's
citations/methods outputs when describing the analysis.
