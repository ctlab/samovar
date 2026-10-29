# Test suites

Run tests using the Python recorded by `samovar.paths.python_path()` after
installation through `./install.sh`:

```bash
python -m pytest -m mandatory
python -m pytest -m optional
python -m pytest --collect-only -q
```

Unmarked tests become `mandatory` in `conftest.py`. Never apply both suite
markers to one test. When a module mixes core behavior and native integrations,
mark the individual integrations rather than making the entire module optional.

## What stays mandatory

- Public prepare/exec wiring, checkpoints, sample/read identities and input validation.
- Annotation/abundance schema, ground truth, merging and splitting paired, SE and hybrid reads.
- Training, persisted model preprocessing, prediction and chunk consistency.
- Tool registration/import, argument forwarding and every interchangeable stage contract.
- Small local numerical and parser checks that can catch a wrong scientific result.

The core suite includes actual CLI, Snakemake, native combiner and tiny ISS
checks. A mock of an optional program checks our adapter contract; it does not
prove that the native program works.

## What stays optional

Native optional programs (R/samovaR, OPAL, MAG toolchains, torch training), live
NCBI/taxdump downloads, larger end-to-end examples, detailed report rendering,
algorithm/resource tuning and tool-specific ecosystem cases. Retain these when
a new tool or a changed adapter could break them. Do not delete a useful
regression merely because the triggering input was unusual.

Delete tests that only pin retired source layouts or copy static constants,
weak presence-only duplicates, and repeated runners that execute the same
baseline contracts already collected directly. Do not replace behavior
assertions with source-text searches.

## Isolation and downloads

Mandatory tests use a private install config and caches. The config keeps the
installed paths of ISS, Snakemake, Kraken2 and Kaiju, but does not inherit the
developer's custom tools, databases or API credentials. Individual tests may
replace that config to test lookup and compatibility behavior.

`tests/data/taxdump` is a deliberately small, NCBI-shaped test taxonomy. It
covers the taxa used by core fixtures, including strain-to-species/genus mapping
and a virus without a genus. It is not a biological reference database.

Python DNS/connection and urllib URL checks reject external requests in
mandatory tests, including requests routed through a local HTTP proxy. Local
HTTP fixtures remain allowed. Child processes inherit the private config and
taxdump; tests launching download-capable external programs must supply local
inputs or be optional. Optional integration tests may use real installed
programs and data; missing dependencies should skip explicitly at execution
time. A skip does not validate that integration.

Use `tmp_path`/temporary directories for outputs. Tests must not delete or
publish into shared `tests_outs/`, `samovar/` run trees or the live install
config. To inspect timing without rerunning the whole suite:

```bash
python -m pytest tests/test_genome_fetcher.py -m mandatory --durations=0
```
