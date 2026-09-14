# Apply a held-out sample, then merge it back

Generate four toy samples, train on three (`samovar exec`), run `samovar apply --full` on the held-out sample, then merge the two runs in both modes and `samovar exec` (checkpoints pick the next stage).

```bash
bash examples/apply_and_merge/pipeline.sh
```

| Directory | What it is |
|-----------|------------|
| `run/generated` | `samovar generate` (4 samples) |
| `run/train` | copy minus 1 sample → prepare (toy Kraken2 + Kaiju) → exec |
| `run/holdout` | copy keeping that 1 sample |
| `run/applied` | `samovar apply --full` of the holdout using models from `train` |
| `run/merged_initial` | `samovar merge --mode initial` then exec (abundance → regen → ML) |
| `run/merged_regenerated` | `samovar merge --mode regenerated` then exec (viz regenerated → reprofile) |

Sample names in the two parents must stay distinct (the split removes the overlapping FASTQ). Downstream exec does not take `--startpoint`: merge writes checkpoints so unfinished stages run.
