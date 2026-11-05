# SAMOVAR generate scenarios

Use inspected genome files and verified taxids. Artifacts belong in tests_outs/.

| Case | Prompt | Expected CLI/refusal and artifacts |
| --- | --- | --- |
| Clear NL | Generate two ISS paired samples with weights 10+10+10 for supplied E. coli (562), yeast (4932) and phage (2886930) FASTAs, 60 pairs/sample, no host, two cores, into RUN. | generate --simulator iss --genomes 562,4932,2886930 --abundance TABLE with those three verified taxids/equal weights; real nonempty paired FASTQs after exec. |
| Taxid list | Use only verified local taxids 562,4932,2886930 with ISS, two samples, 60 pairs, no host and explicit RUN. | Lock the genome list with --genomes, reuse local FASTA, record argv; no accession lookup or download. |
| Ambiguous | Generate 10+10+10. | Ask which organisms/groups and units, simulator and output directory; no generation. |
| Reuse | Re-execute generation for the same known RUN/config and finish its prepared analysis. | Generation checkpoints/read hashes and timestamps preserved. No resimulation/redownload or changed simulator. |
