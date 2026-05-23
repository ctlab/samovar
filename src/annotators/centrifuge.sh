#!/bin/bash
# ==============================================================================
# Script: centrifuge.sh
# Purpose: Snakemake entry point for the Centrifuge classification tool.
# Description: Passes all positional arguments from Snakemake to the custom
#              router script, appending the specific tool flag (-p centrifuge).
# ==============================================================================

bash src/annotators/custom.sh "$@" -p centrifuge
