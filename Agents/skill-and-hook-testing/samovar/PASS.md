# samovar: supported workflow checks passed

Native generated and supplied-read routes both reached viz_reprofiled; each has 120 initial annotation rows, 160 regenerated rows across four samples, a trained model and 120 reprofiled rows. audit.md records data/prompt audit and choices. Ambiguous request refused manually without CLI execution. Partial resume retained completed initial annotation hash/mtime; completed repeat retained all captured files.

Validation date: 2026-11-05. Environment: Linux, configured Python 3.13.2,
SAMOVAR 0.13.1 imported from the working checkout. All five skills were installed
with skills_install.sh into a standalone generic project and read/validated
there. These are manual agent-driven CLI checks in this session, not an
independent model evaluator or separate Cursor/Claude UI sessions.

Evidence root (relative to the SAMOVAR checkout):
`tests_outs/solve-task-122-20261105/`. See `audit.md`, `verification.json`,
`manual-refusals.md`, `negative-cli.json`, the named native command logs and
`generated-ok/` / `supplied/` artifacts. The run config was private. Only
bundled truncated public fixtures and already installed native programs/indexes
were used; no reference download or new example was required.

Actual regeneration uses rescale_abundance=false: each annotator classified
40/60 initial pairs per sample. Its direct tables and eight regenerated mate
files each contain 40 records, not a fabricated 60. Reprofiled predictions cover
all 120 initial records per route. N_reads=60 was the configured budget.

Limitations preserved: #886 (non-1-based abundance sample labels fail in ISS
workflow) and #887 (sort_reads can write through input symlinks). Skills warn
and stop/arrange an owned working copy as appropriate, without patching runtime
code or hiding failed tests. builtin direct emits a missing-BibTeX warning;
that warning is recorded rather than assigning an invented citation.

No lifecycle hooks are installed: the task explicitly permits an instruction-only
MVP. No claim of every simulator, builder, model or IDE integration coverage.
