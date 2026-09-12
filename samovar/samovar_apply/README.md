# Apply a completed pipeline to a new sample

`samovar apply` reuses an already prepared and executed run: same annotators, QC, export/scoring contracts, and trained ML reprofiler. The new input is a **single sample** FASTQ directory.

```bash
bash samovar/samovar_apply/pipeline.sh
```

## Modes

| Command | ML behaviour |
|---------|----------------|
| `samovar apply --input_dir NEW --pipeline SOURCE --output_dir OUT` | Load `SOURCE/reprofiled_annotations/trained_model.joblib` and predict. No fitting. |
| `samovar apply ... --full` | Regenerates abundance tables from the new sample, then **refits** the configured reprofiler (`run_reprofiler`) using the original labeled regenerated annotations. |

`--pipeline` is the original exec directory (must contain `.log/configs/config_init.yaml`). `--output_dir` / `--outdir` / `--directory` is the apply destination; source artifacts are left in place.

This demo builds a tiny dummy-annotator source run (no Kraken/Kaiju databases), then applies it twice.
