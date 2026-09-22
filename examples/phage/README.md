# Phage

Download phage assemblies, `samovar build --index phage_test` when those indexes are missing, `samovar import` them, then `generate` / `prepare` / `exec` twice. Kaiju’s index includes `GCF_000867865.1`; Kraken2’s includes `GCF_000844825.1`. The first run uses `--reindex 1`, the second `--reindex 0`.

```bash
bash examples/phage/pipeline.sh
```

![scores](figures/initial_scores.png)

![F1](figures/initial_F1.png)

[multiqc/multiqc_report.html](multiqc/multiqc_report.html)
