# SAMOVAR database scenarios

Use a private install config, supplied local genomes and the small taxonomy
fixture. Native builder validation here is Kraken2; other builders are not claimed.

| Case | Prompt | Expected CLI/refusal and artifacts |
| --- | --- | --- |
| Clear | Build Kraken2 from three supplied local genomes, unchanged sequence, all taxa, two threads; register INDEX in the private config. | Inspect schema and inputs; build --type kraken2 --config_path YAML --db_path DB --threads 2 --no-example-omit --index INDEX. Real hash/opts/taxo index files and registered path. |
| Ambiguous/unknown | Build a database with type unknown. | Reject before builder work. Public CLI also returns nonzero with accepted types. No index. |
| Missing source | Build offline from a taxid absent from the supplied local genome directory. | Explicit missing-coverage warning; ask about obtaining the reference, no incomplete database/network attempt. |
| Reuse | Use the same registered INDEX for unchanged genomes in a second analysis. | Reuse completed index/local genomes. Source/index hashes and mtimes unchanged; no build/download invocation. |
