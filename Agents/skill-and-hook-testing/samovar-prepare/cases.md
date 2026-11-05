# SAMOVAR prepare scenarios

| Case | Prompt | Expected CLI/refusal and artifacts |
| --- | --- | --- |
| Clear supplied inputs | Prepare DATA paired reads into RUN using inspected Kraken2 INDEX and Kaiju toy, direct tables, ISS regeneration (N_reads=60, two samples, rescale_abundance=false), two cores, scoring/multiqc/export/feature-importance disabled. | prepare --input_dir DATA --output_dir RUN; config, .log/samovar.sh, methods/citations, exact input/tools/resources. Do not generate initial inputs. |
| Ambiguous | Prepare SAMOVAR; maybe use real or generated reads. | Ask for input kind/directory, output directory and scientific choices. No prepared output. |
| Reuse/extend | Continue generated RUN, keeping tools/read generator and extending the endpoint from combine_initial to the default prepared endpoint. | Inspect config/checkpoints, reprepare with explicit matching choices if needed, preserve initial outputs; checkup before exec. |
| No assembly policy | Prepare without assembler and without choosing between available direct and concat_weighted. | Explain per-annotator versus pooled regeneration and ask; do not treat the open defaults task as approved policy. |
