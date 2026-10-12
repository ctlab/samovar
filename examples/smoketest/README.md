# Smoketest

Shortest check that generate, Kraken2, Kaiju, direct regeneration, and reprofiling run.

One sample, 100 reads, no host, no MultiQC. `direct` keeps the observed abundance table instead of drawing a new one. `./install.sh ask` offers this run and asks how many threads to use (default: 3/4 of the CPUs). That count is `samovar prepare --cores` and, when the toy indexes are not already in the install config, `kaiju-mkbwt -n` / `kraken2-build --threads`.

```bash
bash examples/smoketest/pipeline.sh
```

Output: `examples_outdir/smoketest/`.
