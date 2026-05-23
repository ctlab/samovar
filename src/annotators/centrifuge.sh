#!/bin/bash
# Snakemake / CLI entry point: route through custom.sh as -p centrifuge.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
bash "$SCRIPT_DIR/custom.sh" "$@" -p centrifuge
