# Reuse SAMOVAR code and registered tools

**Conditional:** apply only to applied SAMOVAR analyses, including tool/database
selection and run preparation. Development work follows its repository guide.

Use the installed public SAMOVAR CLI and existing stage contracts. Inspect the
installation's tool/database catalog before downloading, building or importing
another copy. Reuse suitable registered objects and existing genome caches.

Reuse the applicable SAMOVAR harness as well: bundled stage skills, shared
rules, installed hooks and public bins. Before writing a helper or recreating a
workflow, check whether that harness already implements the step and invoke it.
Use hooks through their documented entry points when they are installed and
applicable; do not invent hook protocols, enable new lifecycle hooks implicitly,
or duplicate their logic in the agent. Missing reusable CLI/harness behavior is
a reported gap, not permission to build a private replacement pipeline.
For applied runs, call those public commands directly and use their existing
logs or ordinary shell redirection. Do not create a custom command runner just
for logging, retries or stage sequencing when the harness already handles it.
Carry verified catalog/tool selections and shared rules between stage skills;
refresh discovery after an import/build or a changed installation instead of
repeating the same checks for every stage of an unchanged run.

Let SAMOVAR resolve programs and their environments from its main install config.
Do not repair PATH, call a foreign Python, or invent a second config/launcher to
circumvent a failed stage. Use `samovar tools import` and documented `install.sh`
options when registration or installation is actually needed and authorized.

Missing capability, incompatible database or ambiguous catalog selection must
be reported. Ask for the material choice rather than silently replacing a tool,
changing its contract or implementing a substitute analysis in an IDE snippet.
