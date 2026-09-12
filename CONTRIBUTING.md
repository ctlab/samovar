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

Every new code path must have tests covering both:

* its intended logic;
* its compliance with the relevant contracts.

Where a component participates in a pipeline, add integration tests that verify its interaction with the surrounding stages.

Tests belong under `tests/` and must not import implementation code from `examples/`.

For workflow changes, test repeatability and integration where applicable, not only isolated functions.

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
