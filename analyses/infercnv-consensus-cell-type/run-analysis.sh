#!/bin/bash

# This script runs the module across projects
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
scratch_dir="scratch"
ref_dir="references"
mkdir -p ${scratch_dir}
mkdir -p ${ref_dir}

# Create the gene order file for input to inferCNV
Rscript ${script_dir}/00-make-gene-order-file.R \
    --scratch_dir ${scratch_dir} \
    --local_ref_dir ${ref_dir}

# Run individual projects through the module
# Each script prepares the project's normal references and runs inferCNV across relevant project samples

# SCPCP000015: Ewing sarcoma samples
./project-workflows/run-SCPCP000015.sh
