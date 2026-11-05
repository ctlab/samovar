# SAMOVAR execution scenarios

| Case | Prompt | Expected CLI/refusal and artifacts |
| --- | --- | --- |
| Clear | Execute the known prepared RUN to its recorded endpoint. | Inspect config/checkpoints; checkup --output_dir RUN, then exec --output_dir RUN. Verify both annotators and requested stage outputs plus exit status. |
| Missing directory | Resume SAMOVAR. | Ask which prepared run; do not fall through to the CLI current-directory default. |
| Broken tree | Check the supplied empty directory for execution. | checkup --output_dir EMPTY reports missing inputs/nonzero. No exec or fabricated checkpoints. |
| Resume | Continue RUN after combine_initial to the final endpoint without redo, then repeat the completed run. | Initial annotation hash/mtime and completed checkpoints preserved; first sort_reads can change raw byte order. Completed-repeat output hashes/mtimes remain unchanged. |
