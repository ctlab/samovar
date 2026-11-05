# Installer validation — 2026-11-05

Scope: portable bundle installation and the `install.sh` full/ask integration
function. This is **not** a PASS for the SAMOVAR analysis skill or IDE model
selection; no scientific pipeline or real IDE session was exercised here.

Environment: Linux, Bash 5.2.21, Python 3.13.2, pytest 8.4.2.

Executed from the SAMOVAR checkout:

```bash
python -m pytest tests/test_skills_install.py -q
python -m pytest tests/test_install_ask.py -m mandatory -q
bash -n skills_install.sh install.sh
flake8 Agents/install.py tests/test_skills_install.py --count --select=E9,F821,F822,F823 --show-source --statistics
```

Result: **23 passed**, pytest test time 3.23 seconds; shell syntax and critical
lint passed. The standard skill frontmatter validator also accepted `samovar`.
The existing interactive-installer mandatory tests also passed: **28 passed,
5 optional tests deselected** (1.19 seconds). Markdown links and `git diff
--check` passed as well.

| Case | Input / operation | Observed artifacts and behavior |
| --- | --- | --- |
| Happy | Copy only `Agents/` and `skills_install.sh` to a standalone repository; install into projects with spaces in their paths. | Correct Codex, Cursor, Claude and generic skill paths. Shared rules materialized as regular files; still readable after deleting the source. No global instruction files or lifecycle settings created. |
| Ambiguous / conflict | No IDE indicators in explicit mode; foreign same-name skill; modified managed skill. | Nonzero exit and explanatory error. Foreign and edited instructions preserved. `--dry-run` writes nothing. |
| Reuse / update | Repeat installation, then change the source guidance and reinstall an unmodified managed copy. | Repeat yields identical contents; update installs the new source. Unrelated skill preserved. Only manifest-listed skills distributed. |
| Main installer | Invoke the actual `install_agent_harness` function with configured Python and full/ask/plain/opt-out flags. | Full installs into caller's project; non-TTY ask, plain install and opt-out leave it alone. Separate terminal-input tests cover choosing Claude and declining. Full skips CI/unknown IDE. |

Full `install.sh full` was not rerun: this change does not modify package/tool
installation, and reinstalling external bioinformatics environments would not
validate additional harness behavior. Integration evidence above executes its
actual harness function, independently of the heavy installer stages.
