#!/bin/bash
# ==============================================================================
# Script: assembly_hybrid.sh
# Purpose: Snakemake entry point for the Assembly-Hybrid classification tool.
# Description: Passes all positional arguments from Snakemake to the custom
#              router script, appending the specific tool flag (-p assembly_hybrid).
# ==============================================================================

bash src/annotators/custom.sh "$@" -p assembly_hybrid
