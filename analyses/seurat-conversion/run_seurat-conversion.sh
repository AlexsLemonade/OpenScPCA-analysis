#!/bin/bash
set -euo pipefail

# Convert all files in an ScPCA project to Seurat format
# Expects all projects to be found in the `data/current` directory, organized by project id
# Outputs Seurat files to the `results/seurat` directory, organized by project id

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

data_dir="../../data/current"
output_dir="results/seurat"

# get the project as the first arg
project=$1

mkdir -p $output_dir/$project

Rscript scripts/convert-to-seurat.R --input_dir $data_dir/$project --output_dir $output_dir/$project
