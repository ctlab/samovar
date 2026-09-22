Welcome to the SAMOVAR contibution guide! 

Below are listed some guidelines. Please, use pull requests & github issues for the contributions or the  features ideas.

# Contributing to Samovar

Keep changes modular, testable, and consistent with the existing pipeline architecture.

## 1. Contracts and modularity

If a pipeline stage is defined by a clear contract and is intended to be interchangeable:

* define the contract explicitly;
* make implementations conform to it;
* test the contract independently;
* make the pipeline import/integration layer consume the contract rather than a specific implementation;
* keep contracts mutually consistent across pipeline stages.

When introducing a new contract, first check whether existing modules already compose or execute such contracts. Extend or refactor those modules instead of creating parallel mechanisms. `samovar tools import … --pytest` normalizes `--type`, looks it up in `GROUP_TO_TESTNODE`, and runs that one `tests/test_tool_contracts.py::test_*_contract` against `--exec-path` (via `--tool` / `--tool-type`) before writing config; for a new type, document in/out in `CONTRACTS`, map the group, add one `test_*_contract` that calls `_skip_if_other_type` and exercises the dest from `_tool_path`, and ship a baseline that passes it.

## 2. Reuse before duplication

Before implementing new functionality, search for existing modules that already provide the required behavior.

Prefer:

```text
existing module → reuse
existing module → minimal refactor → reuse
```

over:

```text
existing behavior → copy → new implementation
```

Do not duplicate pipeline, configuration, tool, model, or contract logic inside a new command.

## 3. Tests

Tests are split into **mandatory** and **optional** suites (`pytest` markers). The default GitHub Actions workflow (`python-package.yml`) runs **mandatory tests only**:

```text
pytest -m mandatory
```

Optional suite and full suite:

```text
pytest -m optional
pytest
```

**Full integration** (`full-integration.yml`) runs on a published GitHub Release or `workflow_dispatch`. Only `pytest-full` uses `./install.sh full`. Each example job installs the tools that `samovar prepare` records for those pipelines: core (`iss`, `snakemake`) plus `MultiQC`, conda `kraken2`/`kaiju` where the prepared annotator list names them, and `SparseDOSSA2` only on `examples-sparsedossa` (`examples/multiple_tables`, after `libmpfr-dev` and conda `r-rmpfr`). `examples-public` sets `SAMOVAR_CI_LIGHT_INDEXES=1` so realistic/assembly/databases_comparison use locally built `phage_test` indexes instead of downloading standard_8GB/refseq.

Mandatory tests must be fast, deterministic, and cover the core package, contracts, and essential pipeline paths (built-in baselines and small fixtures; no optional programs). Optional tests cover extended integrations, optional dependencies/programs, large datasets, stress/performance, and broader pipeline combinations.

Installation is validated through `install.sh` (the same procedure GitHub Actions uses before pytest). Do not add a second install path for tests.

Test data may come from repository `data/` / `tests/` fixtures or stable public Internet sources. Test code may use `src/`, `tests/`, and installed package APIs. Tests must never import executable code from `examples/`; examples demonstrate the public interface and are not test infrastructure.

New functionality gets mandatory tests if it is core, optional tests if it is ecosystem/extended, or both when both roles apply.

## 4. Repository structure

Follow the existing repository structure and keep responsibilities separated:

* production code belongs in the package;
* tests belong in `tests/`;
* runnable demonstrations belong in `examples/` or the established example locations;
* documentation belongs in the README/Wiki and repository documentation files.

Do not create a new parallel directory or implementation layer when an existing one serves the same purpose.

## 5. Examples and documentation

Examples demonstrate the public `samovar` bin. They are not an implementation layer and not test infrastructure. Read this section before changing `examples/`.

An example is a short sequence of `samovar` commands (`generate`, `prepare`, `exec`, `build`, `import`, `reindex`, `multiqc`, and the other commands in `samovar help`). Do not wrap that sequence in shell functions. Do not call `python -m samovar…` or inline Python for a step the bin already performs.

Catalog names and flags go on the command line. Register a database with `samovar import` (alias of `samovar tools import`). `samovar build --index NAME` records the same catalog entry.

If a step cannot be written as a `samovar` command and needs a large helper, stop. That is missing bin behavior: add the command and a test instead of growing `examples/common.sh` or the example script.

`examples/common.sh` is only for process setup that is not a pipeline step: `PATH`, the output directory, optional `sbatch` around `exec`, and copying figures after the run.

Keep each example small. Changes to public CLI, configuration, contracts, or workflow semantics also update the README or Wiki.

## 6. Validate the whole change

Before considering a change complete:

* run the relevant unit tests;
* run integration tests for affected workflows;
* run relevant existing tests to detect regressions;
* verify that contracts and public interfaces remain consistent;
* update documentation/examples when required.

Prefer extending the existing architecture over introducing special cases for a single workflow.
