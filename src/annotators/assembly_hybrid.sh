#!/bin/bash
# Snakemake / CLI entry point: route through custom.sh as -p assembly_hybrid.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$SCRIPT_DIR/custom.sh" "$@" -p assembly_hybrid
