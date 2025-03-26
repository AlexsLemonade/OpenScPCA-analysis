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
script_dir="scripts"

# Create the gene order file for input to inferCNV
Rscript ${script_dir}/00-make-gene-order-file.R

##### Analysis for SCPCP000015 #####
