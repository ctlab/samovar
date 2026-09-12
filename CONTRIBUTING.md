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

When introducing a new contract, first check whether existing modules already compose or execute such contracts. Extend or refactor those modules instead of creating parallel mechanisms.

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

Tests are split into **mandatory** and **optional** suites (`pytest` markers). GitHub Actions runs **mandatory tests only**:

```text
pytest -m mandatory
```

Optional suite and full suite:

```text
pytest -m optional
pytest
```

Mandatory tests must be fast, deterministic, and cover the core package, contracts, and essential pipeline paths (dummy/small fixtures; no optional programs). Optional tests cover extended integrations, optional dependencies/programs, large datasets, stress/performance, and broader pipeline combinations.

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

Examples demonstrate public Samovar functionality but are not part of the implementation or test infrastructure.

Changes to public behavior, CLI commands, configuration, contracts, or workflow semantics should update the relevant documentation/Wiki.

Keep examples small and focused on demonstrating the intended public interface.

## 6. Validate the whole change

Before considering a change complete:

* run the relevant unit tests;
* run integration tests for affected workflows;
* run relevant existing tests to detect regressions;
* verify that contracts and public interfaces remain consistent;
* update documentation/examples when required.

Prefer extending the existing architecture over introducing special cases for a single workflow.
