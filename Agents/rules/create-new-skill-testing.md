# Validate SAMOVAR skills and hooks

**Conditional:** apply only when adding or changing SAMOVAR skills, rules, hooks
or their installation under `Agents/`.

For each operational skill/hook, keep local cases in the Git-ignored
`Agents/skill-and-hook-testing/<skill-id>/` or CI artifacts: a clear request, an ambiguous request
that must stop with a question/error, and resume/reuse of existing artifacts.
Record each prompt, expected CLI or refusal, required artifacts and notes.

Run the cases against the installed bundle, not only the source Markdown.
Record actual commands, observed outputs, environment and limitations in
`PASS.md` or a CI artifact. Syntax checks and mocked commands do not establish
real SAMOVAR pipeline behavior. Mark unrun cases explicitly; do not fabricate
PASS records. Stage skills are ready only after their relevant real runs pass.
Do not commit these run records or require their presence in the source checkout.

When agent-behavior evaluation is authorized, use fresh subagents with the
installed bundle, natural-language intent and actual input locations. Do not
give them the expected CLI recipe or let them copy an existing example script.
Compare their derived commands and real artifacts against the intended result.
After correcting ambiguous guidance, rerun with a fresh evaluator that has not
seen the diagnosis. Report a stopped or partial trial as such; a corrected
recipe executed by its author is not independent behavioral validation.

Test installer behavior offline in temporary projects: correct discovery paths,
bundled references, repeated installation, preserved user edits, noninteractive
behavior and paths with spaces. No bioinformatics environment install is needed
for installer-only changes. Do not enable an untested lifecycle hook.
