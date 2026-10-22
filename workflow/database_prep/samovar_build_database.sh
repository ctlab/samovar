#!/bin/bash

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
runtime_exports="$("$SCRIPT_DIR/../../bin/samovar" _runtime-env)"
eval "$runtime_exports"

# Function to show usage
show_usage() {
    echo "Usage: samovar build_database [options]"
    echo ""
    echo "Options:"
    echo "  --type TYPE       Type of database to build (kaiju, kraken2, kraken, krakenunique)"
    echo "  --config_path PATH Path to config YAML file"
    echo "  --db_path PATH    Path to store the database"
    echo "  --index NAME      Record this database in the install config"
    echo "  --flags FLAGS     Extra CLI flags stored with --index"
    echo "  --threads N       Workers (kaiju-mkbwt -n, kraken2-build --threads). Same as prepare --cores"
    echo "  --minimizer-len N Kraken/KrakenUniq override (default: 12–15 based on library size)"
    echo "  --example-omit    TEST/TOY ONLY: drop Escherichia from Kraken2 and Phage Phi X from Kaiju"
    echo "  --no-example-omit Keep every taxon even when building from data/test_genomes"
    echo ""
    echo "Example:"
    echo "  samovar build_database --type kaiju --config_path config.yaml --db_path ./database"
}

# Parse arguments
type=""
config_path=""
db_path=""
omit_args=()
index_name=""
db_flags=""
threads="1"
minimizer_args=()

while [[ $# -gt 0 ]]; do
    case $1 in
        --type)
            type="$2"
            shift 2
            ;;
        --config_path)
            config_path="$2"
            shift 2
            ;;
        --db_path)
            db_path="$2"
            shift 2
            ;;
        --example-omit|--no-example-omit)
            omit_args+=("$1")
            shift
            ;;
        --index)
            index_name="$2"
            shift 2
            ;;
        --flags)
            db_flags="$2"
            shift 2
            ;;
        --threads)
            threads="$2"
            shift 2
            ;;
        --minimizer-len)
            minimizer_args=(--minimizer-len "$2")
            shift 2
            ;;
        -h|--help)
            show_usage
            exit 0
            ;;
        *)
            echo "Error: Unknown option '$1'"
            show_usage
            exit 1
            ;;
    esac
done

# Validate required arguments
if [ -z "$type" ] || [ -z "$config_path" ] || [ -z "$db_path" ]; then
    echo "Error: --type, --config_path, and --db_path are required"
    show_usage
    exit 1
fi

# Validate database type
if [ "$type" != "kaiju" ] && [ "$type" != "kraken2" ] && [ "$type" != "kraken" ] && [ "$type" != "krakenunique" ]; then
    echo "Error: --type must be either 'kaiju', 'kraken2', 'kraken', or 'krakenunique'"
    show_usage
    exit 1
fi

# Run snakemake
if ! [[ "$threads" =~ ^[0-9]+$ ]] || [ "$threads" -lt 1 ]; then
    echo "Error: --threads must be a positive integer (got ${threads})"
    exit 1
fi

"$SNAKEMAKE" -s "$SCRIPT_DIR/Snakefile" \
    --configfile "$config_path" \
    --cores "$threads"

# Run the database build script
extra=()
if [ -n "$index_name" ]; then
    extra+=(--index "$index_name")
fi
if [ -n "$db_flags" ]; then
    extra+=(--flags "$db_flags")
fi
"$PYTHON_PATH" "$SCRIPT_DIR/samovar_build_database.py" --type "$type" --config_path "$config_path" --db_path "$db_path" --threads "$threads" "${omit_args[@]}" "${extra[@]}" "${minimizer_args[@]}"
