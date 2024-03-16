#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

data_dir="../../data/current"
output_dir="results/simulated"

project=$1

for dir in $data_dir/$project/*/; do
    sample=$(basename $dir)
    echo "Processing $sample"
    Rscript scripts/simulate-sce.R --sample_dir $dir --output_dir $output_dir/$project/$sample
done

echo "Creating AnnData files"
Rscript scripts/sce-to-anndata.R --dir $output_dir/$project
