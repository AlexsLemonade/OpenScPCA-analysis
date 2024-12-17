#!/bin/bash
set -euo pipefail

# Convert all files in the current data directory to Seurat format
# Outputs Seurat files to the `results/seurat` directory, organized the same as the data directory

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

data_dir="../../data/current"
output_dir="results/seurat"

# in case the output doesn't exist
mkdir -p $output_dir

# run the R script to convert all files
Rscript scripts/convert-to-seurat.R --input_dir $data_dir --output_dir $output_dir
