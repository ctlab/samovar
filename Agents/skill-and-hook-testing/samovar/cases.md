# SAMOVAR umbrella scenarios

Bind DATA, GENOMES, RUN and INDEX to inspected fixtures/registered objects;
record the bindings and argv in the run audit. Never guess a missing binding.

| Case | Prompt | Expected behavior and artifacts |
| --- | --- | --- |
| Generate → analysis | Generate two ISS paired samples from the three supplied genomes, weights 10+10+10, 60 pairs/sample, no host, two cores; use Kraken2 INDEX and Kaiju toy, direct regeneration (N_reads=60, two samples, rescale_abundance=false), then run the supported default endpoint with scoring/reports/export disabled. | Audit, generate, reuse/build only the requested missing index, prepare, checkup, exec. Nonempty initial and regenerated outputs for both annotators, trained reprofiler and final reprofiled outputs. |
| Supplied reads | Run the existing paired reads in DATA with Kraken2 INDEX and Kaiju toy, direct regeneration, same explicit resources and endpoint. Do not simulate initial inputs. | Prepare/checkup/exec without generate, input hashes preserved, both annotators represented and requested later outputs produced. |
| Ambiguous | Generate a metagenome and run SAMOVAR. | Ask for community, simulator, directories and material analysis choices. No generate/prepare/exec and no output directory created. |
| Resume | Continue the known RUN from its completed combine_initial stage to the requested endpoint, then repeat the completed run. | Preserve the completed initial annotation hash/mtime and skip its checkpoints; first sort_reads may change FASTQ byte order. On the already-complete repeat, preserve all captured output hashes/mtimes. |
