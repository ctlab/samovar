#!/bin/bash
# Setup
set -e
cd "$(dirname "$0")/.."
BUILD_DB=false
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -b|--build-db) BUILD_DB=true ;;
        *) echo "[ERROR] Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

if [ -f build/config.json ]; then
    PYTHON_PATH=$(grep -o '"python_path": *"[^"]*"' build/config.json | sed 's/"python_path": *"\(.*\)"/\1/')
    R_PATH=$(grep -o '"r_path": *"[^"]*"' build/config.json | sed 's/"r_path": *"\(.*\)"/\1/')
else
    echo "SamovaR is not installed: check build/config.json"
    exit 1
fi

out_dir="tests_outs"
mkdir -p $out_dir

# ==========================================
# 1. Databases (Conditional Build)
# ==========================================
if [ "$BUILD_DB" = true ]; then
    echo "[INFO] Build flag (-b) detected. Building databases from scratch..."
    snakemake -s workflow/database_prep/Snakefile \
        --configfile workflow/database_prep/config.yaml \
        --cores 16

    $PYTHON_PATH workflow/database_prep/build_database_kraken2.py
    $PYTHON_PATH workflow/database_prep/build_database_kaiju.py
else
    echo "[INFO] Skipping database build. Using pre-built databases from db/ folder."
    echo "       (Run with -b or --build-db to force database generation)"
fi

# ==========================================
# 2. Initial Annotation
# ==========================================
echo "[INFO] Running initial annotators..."
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config_init.yaml \
    --cores 16

echo "[INFO] Combining initial annotation tables..."
$PYTHON_PATH workflow/combine_annotation_tables.py \
    -i tests_outs/benchmarking/initial_reports \
    -o tests_outs/benchmarking/initial_annotations \
    --nodes_dmp db/kraken2/taxonomy/nodes.dmp

echo "[INFO] Visualizing initial annotations..."
$R_PATH -s -f "workflow/compare_annotations.R" \
    --args \
    --annotation_dir tests_outs/benchmarking/initial_annotations \
    --output_dir tests_outs/benchmarking/initial_annotations_plots

# ==========================================
# 3. Read Regeneration
# ==========================================
echo "[INFO] Preparing genomes for regeneration..."
mkdir -p tests_outs/benchmarking/genomes
cp data/test_genomes/meta/* tests_outs/benchmarking/genomes 2>/dev/null || true
cp data/test_genomes/host/* tests_outs/benchmarking/genomes 2>/dev/null || true
cp genomes/mock_community/*.fa tests_outs/benchmarking/genomes 2>/dev/null || true

echo "[INFO] Regenerating synthetic reads..."
snakemake -s workflow/annotation2iss/Snakefile \
    --configfile workflow/annotation2iss/config.yaml \
    --cores 16

echo "[INFO] Cleaning up intermediate fastq files..."
find tests_outs/benchmarking/regenerated -type f -empty -delete 2>/dev/null || true
rm -f tests_outs/benchmarking/regenerated/*_*_*_R*.fastq
rm -f tests_outs/benchmarking/regenerated/*_abundance*
rm -f tests_outs/benchmarking/regenerated/*iss.tmp*

# ==========================================
# 4. Re-Annotation
# ==========================================
echo "[INFO] Running annotators on regenerated reads..."
snakemake -s workflow/annotators/Snakefile \
    --configfile workflow/annotators/config_reannotate.yaml \
    --cores 16

echo "[INFO] Combining regenerated annotation tables..."
$PYTHON_PATH workflow/combine_annotation_tables.py \
    -i tests_outs/benchmarking/regenerated_reports \
    -o tests_outs/benchmarking/regenerated_annotations \
    --nodes_dmp db/kraken2/taxonomy/nodes.dmp

echo "[INFO] Visualizing regenerated annotations..."
$R_PATH -s -f "workflow/compare_annotations.R" \
    --args \
    --annotation_dir tests_outs/benchmarking/regenerated_annotations \
    --output_dir tests_outs/benchmarking/regenerated_annotations_plots

# ==========================================
# 5. ML Feature Extraction
# ==========================================
echo "[INFO] Extracting biological features for ML..."
FEATURES_OUT="tests_outs/benchmarking/features.tsv"
TMP_FASTQ="tests_outs/benchmarking/combined_temporary_R1.fastq"

cat tests_outs/benchmarking/initial/*_R1.fastq > "$TMP_FASTQ"
$PYTHON_PATH src/annotators/fastq_annotator.py "$TMP_FASTQ" -o "$FEATURES_OUT" --chunk_size 50000
rm "$TMP_FASTQ"

# ==========================================
# 6. Machine Learning
# ==========================================
echo "[INFO] Running ML Pipeline..."
$PYTHON_PATH workflow/ML.py \
    --reprofiling_dir tests_outs/benchmarking/initial_annotations \
    --validation_file tests_outs/benchmarking/regenerated_annotations/complex.annotation.csv \
    --features "$FEATURES_OUT" \
    --output_dir tests_outs/benchmarking/ml_results

echo "[INFO] Visualizing final ML results..."
$R_PATH -s -f "workflow/compare_annotations.R" \
    --args \
    --annotation_dir tests_outs/benchmarking/ml_results \
    --output_dir tests_outs/benchmarking/ml_results_plots

echo "[SUCCESS] Pipeline execution complete!"