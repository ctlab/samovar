---
name: samovar-exec
description: Check, execute or resume a known prepared SAMOVAR run, preserving checkpoints and verifying requested outputs rather than relying on exit status alone.
---

# Execute or resume SAMOVAR

Apply only to SAMOVAR run execution, checkup and resume. Read
[scientific integrity](references/rules/scientific-integrity-samovar.md) and
[code reuse](references/rules/code-usage.md).

A directory with `.log/repeats.json` is a package. `samovar exec` on that root runs each `out_XXXX` one after another with the same core and memory budget. A failed repeat does not stop the others. Resume leaves successful checkpoints and seeds alone. In incremental mode, a missing combined annotation table or one missing `*.ingested` sentinel drops every sentinel for that initial stage so annotation runs again for every sample before combine. `--redo` drops those sentinels too. `--redo` applies inside the selected run only. After the children finish, exec rewrites `summary/profiles_long.csv` and `summary/results_long.csv` from saved profiles, scores and runtime. A failed repeat adds no rows. `python -m samovar.repeats tables --root RUN` rebuilds those files without rerunning the children.

Require the run directory from the request or the preceding prepare result.
An unspecified directory is a clarification/error, never the current directory
by assumption. Inspect `samovar checkup --help` and `samovar exec --help` and
the run's configuration, prepared script and `.log/checkpoints/`.

Run `samovar checkup --output_dir RUN` for the intended start/end window before
exec. Inspect the report, not just its exit code. Missing input or a broken tree
must be explained before execution. If the directory contains only generation
configuration, report that state; do not pretend it is a prepared full analysis.

Use the public `samovar exec --output_dir RUN`, adding start/end options only
when requested and supported by the installed help. Default execution resumes
checkpoints. For a paused run, check which stages have completed and which inputs
the next stage needs. Do not delete checkpoints, rewrite provenance, silently
change endpoints or add `--redo` to recover from an error. A redo or changed
analysis needs user authorization; use documented branch facilities when suitable.

Capture the exit status, logs, stage reports and actual requested artifacts.
Verify expected samples/annotators, nonempty outputs where required and any
known-truth/ML/report endpoints. Distinguish a deliberate early endpoint from a
complete analysis and expected empty biology from a failed tool. Report partial
regeneration or missing annotators explicitly even if the command returned zero.

On resume, verify completed stage outputs were reused and requested later stages
were produced. The first `sort_reads` stage can legitimately change FASTQ byte
order/mtime; distinguish that from rerunning a completed generation/annotation.
Check record identities/counts and completed annotation artifacts, and compare
hashes/mtimes on a repeated already-complete run. Return actual paths, coverage, warnings and
unavailable checks, with commands and decisions recorded in the run audit.
