# Stable MAG ids

Toy genomes, a few hundred reads, and the identity assembly slots. After binning, tetranucleotide profiles (k = 4) assign a stable `MAG_ID`. The same catalog is reused when reads are regenerated, so the id is not the bin filename.

```bash
bash examples/mag_identity/pipeline.sh
```

The other feature variant maps reads onto the MAGs saved from the initial annotation and does not keep ids from the new bins:

```bash
SAMOVAR_MAG_MODE=map-initial bash examples/mag_identity/pipeline.sh
```

Both modes use `feat_*_MAG_ID`, with the same column on initial and regenerated annotations. The `map-initial` mode requires minimap2 and samtools (`./install.sh minimap2`); the identity taxonomy in this example is a debugging baseline, not a GTDB gene check. Production GTDB mapping uses `--mag-taxonomy gtdbtk` and a configured GTDB database.
