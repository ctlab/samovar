# Generate reads with InSilicoSeq
snakemake -s workflow/iss_test/Snakefile \
    --configfile workflow/iss_test/config.yaml \
    --cores 1

# optional: build custom databases
if true; then
    # Subset genomes for database creation
    snakemake -s workflow/database_prep/Snakefile \
        --configfile workflow/database_prep/config.yaml \
        --cores 1

    # Prepare databases
    python workflow/database_prep/build_database_kraken2.py
    python workflow/database_prep/build_database_kaiju.py
fi

# Run annotators on initial reads
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config_init.yaml \
    --cores 1

# Combine annotation tables
python workflow/combine_annotation_tables.py \
    -i tests_outs/benchmarking/initial_reports \
    -o tests_outs/benchmarking/initial_annotations

# Visualize annotations
Rscript workflow/compare_annotations.R \
    --annotation_dir tests_outs/benchmarking/initial_annotations \
    --output_dir tests_outs/benchmarking/initial_annotations_plots

# Add pre-downloaded genomes to the genome directory
mkdir -p tests_outs/benchmarking/genomes
cp data/test_genomes/meta/* tests_outs/benchmarking/genomes
cp data/test_genomes/host/* tests_outs/benchmarking/genomes

# Translate annotation table to new reads set
snakemake -s workflow/annotation2iss/Snakefile \
    --configfile workflow/annotation2iss/config.yaml \
    --cores 1

# Clean empty annotation files
find tests_outs/benchmarking/regenerated -type f -empty -delete

# Run annotators on new reads set
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config_reannotate.yaml \
    --cores 1

# Combine annotation tables
python workflow/combine_annotation_tables.py \
    -i tests_outs/benchmarking/regenerated_reports \
    -o tests_outs/benchmarking/regenerated_annotations

# Visualize & combine results
Rscript workflow/compare_annotations.R \
    --annotation_dir tests_outs/benchmarking/regenerated_annotations \
    --output_dir tests_outs/benchmarking/regenerated_annotations_plots \
    --csv tests_outs/benchmarking/regenerated_annotations/combined_annotation_table.csv
