# Operational cases — not yet run

These cases are for a real agent session with the installed skill. Replace
fixture placeholders with inspected paths before execution. A transcript must
show actual CLI and artifact checks; installer pytest is separate evidence.

| Case | Prompt / inputs | Expected behavior | Artifacts / notes |
| --- | --- | --- | --- |
| Clear | Annotate this supplied tiny paired-read dataset with these two registered toy databases; prepare and execute into this fresh output directory using two cores. | Inspect files/catalog/help, prepare then checkup/exec using supported flags; verify both annotators and requested final stages. | Commands, audit, prepared configs, checkpoints, nonempty results for both tools. Requires real tools and known fixtures. |
| Ambiguous | Generate a metagenome and run SAMOVAR. | Ask for community, simulator, output location and material analysis choices; do not start generate or prepare. | Transcript with clarification and no generated data. |
| Reuse/resume | Continue this supplied prepared, partially completed tiny run with its existing registered databases. | Inspect config/checkpoints/checkup, resume the same directory; reuse completed work without redownload or blanket redo. | Before/after checkpoint and output evidence; real resume transcript. |
