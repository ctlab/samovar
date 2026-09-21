# All annotators

Ten NCBI bacterial genomes and four classifiers: Kraken2 `standard_8GB`, Kaiju `refseq`, KrakenUniq `microbial`, and Kraken 1 `minikraken_4GB`.

Each name is a catalog entry with an official download URL. If that index (or another installed version of the same tool) is already in the SamovaR database catalog, the run uses it. Otherwise the URL is downloaded under `examples_outdir/databases/`.

```bash
bash examples/all_annotators/pipeline.sh
```

Large indexes: `SAMOVAR_SLURM=1 SAMOVAR_SLURM_CPUS=16 SAMOVAR_SLURM_MEM=480G`.

Figures land in `figures/`, and the MultiQC report in `multiqc/multiqc_report.html`. The full run tree is `examples_outdir/all_annotators/`.
