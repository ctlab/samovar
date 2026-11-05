# Smoketest

Shortest check that generate, Kraken2, Kaiju, direct regeneration, and reprofiling run.

One sample, 100 reads, two genomes (4932 and 562), no host, no MultiQC. `concat` merges Kraken2 and Kaiju into one abundance table, so regeneration writes one FASTQ rather than one sample per annotator. ISS stays on one CPU: at this depth the error-model load is single-threaded, and more workers only add empty temp files. `./install.sh ask` offers this run and asks how many threads to use (default: 3/4 of the CPUs). That count is `samovar prepare --cores` and, when the toy indexes are not already in the install config, `kaiju-mkbwt -n` / `kraken2-build --threads`.

```bash
bash examples/smoketest/pipeline.sh
```

Output: `examples_outdir/smoketest/`.
