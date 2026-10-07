#!/bin/bash
# Snakemake / CLI entry point: route through custom.sh as -p metax.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$SCRIPT_DIR/custom.sh" "$@" -p metax
