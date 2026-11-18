---
name: samovar-tools-import
description: Turn a plain-language description of an external program into a SAMOVAR tools import, with contract questions, a tool-specific pytest, and a smoke check in the target stage.
---

# Import a program into SAMOVAR

Use this skill when someone describes a program that should become a
`samovar tools import` tool. The public command is `samovar tools import`
(`samovar import` is the same command). `--pytest` runs the shipped contract
node for that group. It does not write a new test.

Do not modify `src/` or `tests/` of the SAMOVAR checkout to make an import
pass, and do not weaken a contract. A missing `--type` or contract test is a
gap report, not a new group created from this skill. Do not rewrite built-in
reprofilers (`rfc2`, `random_forest`, `ensemble`, `adaboost`, `linear`).

## 1. Read the description

Take the natural-language description, and any path, `--help`, example
command, or guessed stage the user already gave. Read `samovar tools --status`
and `samovar tools import --help` before choosing a group.

Current groups are the `TOOL_GROUPS` list (`annotator`, `profiler`,
`reads_generator`, `metagenome_generator`, `table_reads_generator`,
`table_scoring`, `concordance`, `sample_scoring`, `sample_filtering`,
`scoring`, `reprofiler`, `annotation_converter`, `export`, `qc`, `assembler`,
`gene_caller`, `binner`, `binner_qc`, `binner_combine`, `mag_taxonomy`,
`aligner`, `mag_quantifier`, `taxon_quantifier`, `read_assigner`,
`mag_homology`, `mag_mapping`, `feature_importance`). `runtime`, `compiler`,
`workflow`, and `database` are not analysis stages for this skill.

If the description fits no group, or fits two groups with different pipeline
meanings, stop and say so. Do not pick one to make the import succeed.

## 2. Ask only what changes the contract

Ask only what the description left open and what would change the import or
the scientific result:

- what the program reads and writes
- per-read labels, a sample profile, ML correction, trimming, table
  regeneration, or scoring
- whether a database is required and whether the binary is installed
- for ML, whether a Python `reprofile()` exists or the program is CLI-only
- a name clash with a tool already in `samovar tools --status`
- a BibTeX citation, when the user wants one recorded

Do not ask for facts already stated. `python3` (or `python`) offered as the
reprofiler itself is a refusal: the interpreter is not a model. Say that and
stop.

## 3. Choose native mapping or an adapter

Prefer the native group when SAMOVAR already launches that kind of program
with its own flags, including QC (`--type qc`) and annotators or profilers
that already have an adapter. Use `--flags-translate` when the binary's own
flags must be mapped onto SAMOVAR's `-i`/`-o`/`--mem` contract. Do not invent
a Python wrapper for a native tool such as fastp.

Use an adapter plus a manifest when the program does not speak the stage
contract. The adapter is a Python file in the user's project, with a
documented entry point. The manifest records `name`, `type`, `exec-path`,
`flags`, `flags-translate`, and `bibtex` when present. Show the exact
`samovar tools import` command as a dry run. Do not run it yet.

A CLI-only ML program that has no `reprofile()` will fail `--pytest` for
`--type ml` / `reprofiler`. Say that, and offer a Python wrapper that
implements `reprofile()`. Import without `--pytest` only after the user
explicitly accepts it, and label the result "not contract-certified".

## 4. Tests, then the stage

Create the tool-specific pytest beside the adapter in the user's project, not
under SAMOVAR `tests/`. On a small fixture it must call the real program or
adapter and check non-empty, meaningful output (columns, taxids, row or read
counts), not only the exit code. Include negative cases: empty or corrupt
input, wrong columns, and a missing database when the tool needs one. Those
must fail with the expected error, not an empty success.

Then run the shipped contract:

```bash
samovar tools import -n NAME -t TYPE --exec-path PATH --pytest
```

The node comes from `GROUP_TO_TESTNODE` (for example
`tests/test_tool_contracts.py::test_reprofiler_contract`). Use a temporary
`SAMOVAR_CONFIG` for this check.

The stage smoke is required, not optional. Use the same function the pipeline
uses:

| Group | Smoke |
|---|---|
| `reprofiler` | `prepare --reprofiler NAME`, then the reprofile stage. The table is a non-empty `*_reprofiled.csv` with `taxid_SAMOVAR`. |
| `annotator` | A TSV with `seq` plus a tax or feature column. |
| `profiler` | `taxid`, `N_*`, and `.profile.json`. |
| `qc` | A non-empty destination FASTQ. |
| other | The stage function that consumes that group. |

`fastp` is `--type qc` with native mapping and no fake Python wrapper. A
logistic classifier described as an ML program is `--type ml` and must provide
`reprofile()` before it is called contract-certified.

## 5. Write the config only after yes

Reading `--help` and running pytest against a temporary config do not need
confirmation. Writing the install config and installing packages do. Run the
import command only after the user confirms that exact command.

If the work is blocked (no executable, wrong group, contract failed, native
path unsupported, or the user said no), report that and stop. Do not swap in
another tool or change `--type` so the import succeeds.

## Session files

Keep the session in the user's project so it can be repeated with
`bash commands.sh` on the same SAMOVAR version: `description.md`,
`questions.md`, `manifest.yaml`, `adapter.py` when needed, `test_<tool>.py`,
`fixtures/`, `commands.sh`, `pytest.log`, `contract_pytest.log`, `smoke.log`,
and `approval.txt`. Do not put secrets or machine-local trees in those files.
