# Validation records

The installer is tested by `tests/test_skills_install.py` using temporary
projects and no network or installed bioinformatics tools. That suite also gates
the presence of scenarios and evidence for each manifest-listed skill. Operational
cases and their actual native/manual validation are recorded in the five
per-skill directories. See each `PASS.md` for coverage and limitations; filesystem
installation tests do not establish scientific correctness or IDE UI behavior.

Follow [the testing rule](../rules/create-new-skill-testing.md). Keep large run
artifacts under `tests_outs/`, and only concise prompts, expectations and real
evidence summaries here. Never label a planned case PASS.
