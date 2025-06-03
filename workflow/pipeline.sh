# Setup
set -e

if [ -f build/config.json ]; then
    PYTHON_PATH=$(grep -o '"python_path": *"[^"]*"' build/config.json | sed 's/"python_path": *"\(.*\)"/\1/')
    R_PATH=$(grep -o '"r_path": *"[^"]*"' build/config.json | sed 's/"r_path": *"\(.*\)"/\1/')
    R_LIB_PATH=$(grep -o '"r_lib_path": *"[^"]*"' build/config.json | sed 's/"r_lib_path": *"\(.*\)"/\1/')
else
    echo "SamovaR is not installed: check build/config.json"
    exit 1
fi

out_dir="tests_outs"
mkdir -p $out_dir

# optional: build custom databases
if true; then
    # Subset genomes for database creation
    snakemake -s workflow/database_prep/Snakefile \
        --configfile workflow/database_prep/config.yaml \
        --cores 1

    # Prepare databases
    $PYTHON_PATH workflow/database_prep/build_database_kraken2.py
    $PYTHON_PATH workflow/database_prep/build_database_kaiju.py
fi

# optional: generate reads with InSilicoSeq for automated benchmarking;
# otherwise, use real data
if true; then
    snakemake -s workflow/iss_test/Snakefile \
        --configfile workflow/iss_test/config.yaml \
        --cores 1
fi

# Run annotators on initial reads
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config_init.yaml \
    --cores 1

# Combine annotation tables
$PYTHON_PATH workflow/combine_annotation_tables.py \
    -i tests_outs/benchmarking/initial_reports \
    -o tests_outs/benchmarking/initial_annotations

# Visualize annotations
$R_PATH -s -f "workflow/compare_annotations.R" \
    --args \
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

# Clean up
find tests_outs/benchmarking/regenerated -type f -empty -delete
rm tests_outs/benchmarking/regenerated/*_*_*_R*.fastq
rm tests_outs/benchmarking/regenerated/*_abundance*
rm tests_outs/benchmarking/regenerated/*iss.tmp*

# Run annotators on new reads set
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config_reannotate.yaml \
    --cores 1

# Combine annotation tables
$PYTHON_PATH workflow/combine_annotation_tables.py \
    -i tests_outs/benchmarking/regenerated_reports \
    -o tests_outs/benchmarking/regenerated_annotations \
    -s 2

# Visualize & combine results
$R_PATH -s -f "workflow/compare_annotations.R" \
    --args \
    --annotation_dir tests_outs/benchmarking/regenerated_annotations \
    --output_dir tests_outs/benchmarking/regenerated_annotations_plots \
    --csv tests_outs/benchmarking/regenerated_annotations/combined_annotation_table.csv

# Train and test ML
$PYTHON_PATH workflow/ML.py \
    --reprofiling_dir tests_outs/benchmarking/initial_annotations \
    --validation_file tests_outs/benchmarking/regenerated_annotations/combined_annotation_table.csv \
    --output_dir tests_outs/benchmarking/reprofiled_annotations

# Check reprofiled results
$R_PATH -s -f "workflow/compare_annotations.R" \
    --args \
    --annotation_dir tests_outs/benchmarking/reprofiled_annotations \
    --output_dir tests_outs/benchmarking/reprofiled_annotations_plots \
    --csv tests_outs/benchmarking/reprofiled_annotations/combined_annotation_table.csv