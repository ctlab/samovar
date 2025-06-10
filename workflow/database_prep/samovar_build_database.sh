#!/bin/bash

set -e

# Function to show usage
show_usage() {
    echo "Usage: samovar build_database [options]"
    echo ""
    echo "Options:"
    echo "  --type TYPE       Type of database to build (kaiju, kraken2, kraken, krakenunique)"
    echo "  --config_path PATH Path to config YAML file"
    echo "  --db_path PATH    Path to store the database"
    echo ""
    echo "Example:"
    echo "  samovar build_database --type kaiju --config_path config.yaml --db_path ./database"
}

# Parse arguments
type=""
config_path=""
db_path=""

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
snakemake -s $(dirname "$0")/Snakefile \
    --configfile $config_path \
    --cores 1

# Run the database build script
python "$(dirname "$0")/samovar_build_database.py" --type "$type" --config_path "$config_path" --db_path "$db_path" 