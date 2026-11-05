# Portable SAMOVAR agent bundle

This is the **distribution source** for new SAMOVAR analysis skills. The existing
ignored `.agents/` directory serves development of this checkout and is never
copied into this bundle. No host-specific Vikunja client or old skills are shipped.

## Install from an IDE terminal

From the project where you want the skills (Python 3.10+ and Bash):

```bash
bash /path/to/samovar/skills_install.sh
```

An active Codex, Cursor Agent or Claude Code session is detected from its
environment. An ordinary IDE terminal may not expose its agent; choose explicitly:

```bash
bash /path/to/samovar/skills_install.sh --ide cursor
bash /path/to/samovar/skills_install.sh --ide codex --project /path/to/analysis
bash /path/to/samovar/skills_install.sh --ide claude --dry-run
```

| Target | Project-local installation |
| --- | --- |
| Codex | `.agents/skills/<name>/` |
| Cursor | `.cursor/skills/<name>/` |
| Claude Code | `.claude/skills/<name>/` |
| Other IDE (`--ide generic`) | `samovar-harness/skills/<name>/`; point the agent at `SKILL.md` or copy to its documented skill directory |

Reload the IDE's skill list/session if needed, then ask it to use `samovar` for
your analysis. Generic installation is a portable bundle, not a claim that every
IDE automatically discovers the same directory. Discovery paths follow the
[Codex](https://learn.chatgpt.com/docs/build-skills),
[Cursor](https://cursor.com/docs/skills) and
[Claude Code](https://code.claude.com/docs/en/skills) documentation.

`install.sh full` attempts this installation when an IDE session is detected;
it skips CI and unknown terminals. `install.sh ask` (including `full ask`) asks
which target to install; noninteractive ask skips it. The destination is the
directory from which `install.sh` was invoked. Set `SAMOVAR_HARNESS_PROJECT` to
override it, `SAMOVAR_IDE` to select a target, or `SAMOVAR_INSTALL_HARNESS=0` to
opt out. Plain `install.sh` and tool-only installations do not install skills.
The harness uses the same configured Python when called by `install.sh`.

Re-running the command updates unmodified managed skills. A same-name foreign
skill or local edits cause an error before changes; preserve/rename that copy
before retrying. Other skills, settings and global instruction files are left
alone. No lifecycle hooks are enabled. Installation does not download anything
or invoke a bioinformatics pipeline.

## Publishing boundary and layout

The portable payload is **`Agents/` plus root `skills_install.sh`**. These paths
can be copied unchanged into a standalone skills repository; the installer has
no dependency on `src/samovar`, the package config, or `.agents/`. Publish only
after review; this change does not push or commit the main repository.

- `bundle.json`: explicit list of installable skills: `samovar`,
  `samovar-generate`, `samovar-build-database`, `samovar-prepare`, `samovar-exec`.
- `skills/<name>/SKILL.md`: standard skill format, without IDE-specific APIs.
- `rules/*.md`: canonical conditional SAMOVAR policies. The source skills'
  `references/rules` link points here. Installation materializes regular files
  inside each skill so it is self-contained and survives moving/deleting the
  source checkout. New skills needing rules use the same relative link.
- `hooks/`: hook policy; no event handlers or lifecycle configuration yet.
- `skill-and-hook-testing/`: cases and honest operational validation records.
- `install.py`: shared installer implementation; the root script only launches it.

Rules activate through relevant skills, not unconditional global instructions.
There is no shared native conditional-rule format across IDEs. Bundled Markdown
with an explicit SAMOVAR condition is the portable contract. Do not install
these files into an IDE's always-loaded rules directory. Any future native
adapter must preserve the same condition (e.g. Cursor `alwaysApply: false`).

## Operational scope and validation

`samovar` audits the inputs/request and routes to the four stage skills. They
use the installed CLI/help and registered tools rather than embedding a second
pipeline. Sequential execution works in any agent; optional delegation requires
support and user authorization. Hooks in the Skills tasks are optional; this
bundle ships the instruction-only implementation.

Per-skill scenarios and evidence live in `skill-and-hook-testing/<name>/`.
The mandatory packaging check requires scenarios and a `PASS.md` for every
manifest-listed skill, so a new skill without records fails CI. That presence
check is not a scientific validator: reviewers must inspect the recorded real
commands, outputs and limitations. Installer tests, manual agent checks and
native pipeline execution are reported separately; no claim is made that every
model or IDE UI has been tested.

The current ISS CLI has a known sample-label limitation for supplied abundance
tables (follow-up #886). The generate skill explains when to stop and warn;
it does not silently rename samples or treat a failed run as successful.

When extending this bundle, read the conditional
[layout](rules/agents-portable-layout.md) and
[testing](rules/create-new-skill-testing.md) rules. Operational skills must also
reference [scientific integrity](rules/scientific-integrity-samovar.md) and
[code reuse](rules/code-usage.md). Add the skill name to `bundle.json` only when
it is intended for distribution. Do not copy current checkout-maintenance skills.

Validate installer changes with `python -m pytest tests/test_skills_install.py`
in the SAMOVAR checkout. For a standalone copy use `bash skills_install.sh
--ide generic --project /path/to/temporary/project --dry-run`, then perform a
real install there and inspect its bundled references. This does not replace
the per-skill scientific workflow checks.
