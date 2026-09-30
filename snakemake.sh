#!/bin/bash
# Local run: ./snakemake.sh [snakemake args], e.g. ./snakemake.sh -n for a dry run.
set -eo pipefail

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate snake

snakemake --workflow-profile workflow/profiles/default -s workflow/Snakefile "$@"
