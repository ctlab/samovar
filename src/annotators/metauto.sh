#!/bin/bash
# ==============================================================================
# Script: metauto.sh
# Purpose: Snakemake entry point for the Metauto-DL classification tool.
# Description: Passes all positional arguments from Snakemake to the custom
#              router script, appending the specific tool flag (-p metauto).
# ==============================================================================

bash src/annotators/custom.sh "$@" -p metauto
