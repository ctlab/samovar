# Generate reads with InSilicoSeq
snakemake -s workflow/iss_test/Snakefile \
    --configfile workflow/iss_test/config.yaml \
    --cores 1

# Subset genomes for database creation
snakemake -s workflow/database_prep/Snakefile \
    --configfile workflow/database_prep/config.yaml \
    --cores 1

# Prepare databases
python workflow/database_prep/build_database_kraken2.py
python workflow/database_prep/build_database_kaiju.py

# Run annotators on initial reads
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config_init.yaml \
    --cores 1

# Combine annotation tables
python workflow/combine_annotation_tables.py \
    -i tests_outs/benchmarking/initial_reports \
    -o tests_outs/benchmarking/initial_annotations \
    -t "[^-]*"

# Translate annotation table to new reads set
snakemake -s workflow/annotation2iss/Snakefile \
    --configfile workflow/annotation2iss/config.yaml \
    --cores 1

# Run annotators on new reads set
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config_reannotate.yaml \
    --cores 1

# Combine annotation tables
python workflow/combine_annotation_tables.py \
    -i tests_outs/benchmarking/regenerated_reports \
    -o tests_outs/benchmarking/regenerated_annotations \
    -t "[^-]*" -s 2