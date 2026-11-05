# Reuse SAMOVAR code and registered tools

**Conditional:** apply only to applied SAMOVAR analyses, including tool/database
selection and run preparation. Development work follows its repository guide.

Use the installed public SAMOVAR CLI and existing stage contracts. Inspect the
installation's tool/database catalog before downloading, building or importing
another copy. Reuse suitable registered objects and existing genome caches.

Let SAMOVAR resolve programs and their environments from its main install config.
Do not repair PATH, call a foreign Python, or invent a second config/launcher to
circumvent a failed stage. Use `samovar tools import` and documented `install.sh`
options when registration or installation is actually needed and authorized.

Missing capability, incompatible database or ambiguous catalog selection must
be reported. Ask for the material choice rather than silently replacing a tool,
changing its contract or implementing a substitute analysis in an IDE snippet.
