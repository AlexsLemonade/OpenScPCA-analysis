#!/bin/bash

# This script runs the module workflow.
#
# Usage:
#
# ./run-analysis.sh

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define directories
notebook_dir="notebooks"

##### Analysis for SCPCP000015 #####

# Explore distribution of cell types
Rscript -e "rmarkdown::render('${notebook_dir}/01_ewings-consensus-cell-types.Rmd')"
