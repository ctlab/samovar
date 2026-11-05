# Validation records

The installer is tested by `tests/test_skills_install.py` using temporary
projects and no network or installed bioinformatics tools. The `samovar/` cases
describe the next operational checks; they have **not** been executed as agent
end-to-end tests. There is no operational PASS claim yet.

Follow [the testing rule](../rules/create-new-skill-testing.md). Keep large run
artifacts under `tests_outs/`, and only concise prompts, expectations and real
evidence summaries here. Never label a planned case PASS.
