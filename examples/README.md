# Examples

Each folder is one pipeline (`pipeline.sh` plus a short README). From any working directory:

```bash
bash examples/toy/pipeline.sh
```

Large jobs can use the cluster without changing the scripts: `SAMOVAR_SLURM=1 SAMOVAR_SLURM_CPUS=16`.

| Example | What it shows |
|---------|----------------|
| [toy](toy/) | Basic generate → prepare → exec |
| [scoring](scoring/) | Choosing and comparing annotators |
| [reprofiling](reprofiling/) | A custom ML reprofiler |
| [logistic-correction](logistic-correction/) | Single-annotator logistic abundance correction |
| [multiple_tables](multiple_tables/) | Scoring several abundance-table methods |
| [phage](phage/) | Named SamovaR databases and `reindex` |
| [databases_comparison](databases_comparison/) | Same community, several Kraken2 indexes |
| [realistic](realistic/) | NCBI genomes with public Kraken2 and Kaiju indexes |

Figures and MultiQC reports are copied next to each `pipeline.sh` after a run.
