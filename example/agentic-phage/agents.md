# Agentic phage benchmark

`/samovar` runs this profile using the installed SAMOVAR harness. Paths are
relative to the containing checkout; honor an operator's output override.

Compare two phage communities using ISS paired reads, Kraken2 and Kaiju, then
regenerate communities, train the `random_forest` reprofiler and reprofile both initial
communities. Produce scores and a MultiQC report for each run.

| Input | Assembly accessions |
| --- | --- |
| Community 1 | GCF_000819615.1, GCF_000840245.1, GCF_000836945.1, GCF_000844825.1 |
| Community 2 | GCF_000819615.1, GCA_035764635.1, GCF_000836945.1, GCF_000844825.1, GCF_000867865.1 |
| Kaiju reference | GCF_000819615.1, GCF_000840245.1, GCF_000836945.1, GCF_000867865.1 |
| Kraken2 reference | GCF_000840245.1, GCF_000836945.1, GCF_000844825.1 |

Use `data/test_genomes/host/9606.fna` as 15% host. Each community has two samples
of 2,000 paired fragments, seed 42 and two cores. Regenerate two samples using
the same nominal read budget and default mass-preserving `concat_weighted`;
actual regenerated depths may differ and must be reported. Keep the installed
pipeline defaults for other settings. Do not add assembly or substitute tools.
For an explicitly requested smoke test, use 200 instead of 2,000 for initial
reads and regenerated depth; everything else stays the same.

Reuse matching `phage_test` indexes and cached genomes; downloading these public
assemblies and building missing references through SAMOVAR is authorized.
Register the indexes as `phage_test`. Index first-community downloads in the
shared catalog; keep second-community raw downloads with its run. Retain normal
pipeline catalog harvesting, which can index them later. Ask before installing
missing software.

Save runs under `examples_outdir/agentic-phage/database1` and `database2`, with
reference build data under `examples_outdir/agentic-phage/database`. Resume
matching work and report verified artifacts and warnings.
