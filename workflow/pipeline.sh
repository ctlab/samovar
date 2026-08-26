# Demo pipeline (repo-relative). Prefer `samovar prepare` / `samovar exec` for real runs.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export SAMOVAR_ROOT="$ROOT"
export PATH="$ROOT/bin:$PATH"

# Load install env from build/config_path (custom config) or default XDG.
XDG_CONFIG_HOME="${XDG_CONFIG_HOME:-$HOME/.config}"
_SAMOVAR_CFG=""
if [ -f "$ROOT/build/config_path" ]; then
    _SAMOVAR_CFG="$(tr -d '[:space:]' < "$ROOT/build/config_path")"
fi
if [ -n "$_SAMOVAR_CFG" ]; then
    export SAMOVAR_CONFIG="$_SAMOVAR_CFG"
    if [ -f "$(dirname "$_SAMOVAR_CFG")/env" ]; then
        # shellcheck disable=SC1090
        . "$(dirname "$_SAMOVAR_CFG")/env"
    fi
elif [ -f "$XDG_CONFIG_HOME/samovar/env" ]; then
    # shellcheck disable=SC1090
    . "$XDG_CONFIG_HOME/samovar/env"
fi

CFG_JSON="${SAMOVAR_CONFIG:-}"
if [ -z "$CFG_JSON" ] || [ ! -f "$CFG_JSON" ]; then
    if [ -n "$_SAMOVAR_CFG" ] && [ -f "$_SAMOVAR_CFG" ]; then
        CFG_JSON="$_SAMOVAR_CFG"
    elif [ -f "$XDG_CONFIG_HOME/samovar/config.json" ]; then
        CFG_JSON="$XDG_CONFIG_HOME/samovar/config.json"
    elif [ -f "$ROOT/build/config.json" ]; then
        CFG_JSON="$ROOT/build/config.json"
    fi
fi
unset _SAMOVAR_CFG

if [ -n "${CFG_JSON:-}" ] && [ -f "$CFG_JSON" ]; then
    PYTHON_PATH=$(python3 -c "
import json,sys
p=json.load(open(sys.argv[1]))
print((p.get('compilers') or {}).get('python') or p.get('python_path') or '')
" "$CFG_JSON" 2>/dev/null || true)
fi
PYTHON_PATH=${PYTHON_PATH:-python3}

out_dir="tests_outs"
mkdir -p $out_dir

# optional: build custom databases
if true; then
    snakemake -s "$ROOT/workflow/database_prep/Snakefile" \
        --configfile "$ROOT/workflow/database_prep/config.yaml" \
        --cores 1

    $PYTHON_PATH "$ROOT/workflow/database_prep/build_database_kraken2.py" \
        --config_path "$ROOT/workflow/database_prep/config.yaml"
    $PYTHON_PATH "$ROOT/workflow/database_prep/build_database_kaiju.py" \
        --config_path "$ROOT/workflow/database_prep/config.yaml"
fi

if true; then
    snakemake -s "$ROOT/workflow/iss_test/Snakefile" \
        --configfile "$ROOT/workflow/iss_test/config.yaml" \
        --cores 1
fi

snakemake -s "$ROOT/workflow/annotators/Snakefile" \
    --configfile "$ROOT/workflow/annotators/config_init.yaml" \
    --cores 1

$PYTHON_PATH "$ROOT/workflow/combine_annotation_tables.py" \
    -i tests_outs/benchmarking/initial_reports \
    -o tests_outs/benchmarking/initial_annotations

$PYTHON_PATH "$ROOT/workflow/compare_annotations.py" \
    --annotation_dir tests_outs/benchmarking/initial_annotations \
    --output_dir tests_outs/benchmarking/initial_annotations_plots

mkdir -p tests_outs/benchmarking/genomes
TEST_GENOMES="${SAMOVAR_TEST_GENOMES:-$ROOT/data/test_genomes}"
if [ -d "$TEST_GENOMES/meta" ]; then
    cp "$TEST_GENOMES/meta/"* tests_outs/benchmarking/genomes/ 2>/dev/null || true
fi
if [ -d "$TEST_GENOMES/host" ]; then
    cp "$TEST_GENOMES/host/"* tests_outs/benchmarking/genomes/ 2>/dev/null || true
fi

snakemake -s "$ROOT/workflow/annotation2iss/Snakefile" \
    --configfile "$ROOT/workflow/annotation2iss/config.yaml" \
    --cores 1

find tests_outs/benchmarking/regenerated -type f -empty -delete || true
rm -f tests_outs/benchmarking/regenerated/*_*_*_R*.fastq || true
rm -f tests_outs/benchmarking/regenerated/*_abundance* || true
rm -f tests_outs/benchmarking/regenerated/*iss.tmp* || true

snakemake -s "$ROOT/workflow/annotators/Snakefile" \
    --configfile "$ROOT/workflow/annotators/config_reannotate.yaml" \
    --cores 1

$PYTHON_PATH "$ROOT/workflow/combine_annotation_tables.py" \
    -i tests_outs/benchmarking/regenerated_reports \
    -o tests_outs/benchmarking/regenerated_annotations \
    -s 2

$PYTHON_PATH "$ROOT/workflow/compare_annotations.py" \
    --annotation_dir tests_outs/benchmarking/regenerated_annotations \
    --output_dir tests_outs/benchmarking/regenerated_annotations_plots \
    --csv tests_outs/benchmarking/regenerated_annotations/combined_annotation_table.csv

$PYTHON_PATH "$ROOT/workflow/ML.py" \
    --reprofiling_dir tests_outs/benchmarking/initial_annotations \
    --validation_file tests_outs/benchmarking/regenerated_annotations/combined_annotation_table.csv \
    --output_dir tests_outs/benchmarking/reprofiled_annotations

$PYTHON_PATH "$ROOT/workflow/compare_annotations.py" \
    --annotation_dir tests_outs/benchmarking/reprofiled_annotations \
    --output_dir tests_outs/benchmarking/reprofiled_annotations_plots \
    --csv tests_outs/benchmarking/reprofiled_annotations/combined_annotation_table.csv
