#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

data_dir="../../data/current"
output_dir="results/seurat"

project=$1

mkdir -p $output_dir/$project

Rscript scripts/convert-seurat.R --dir $data_dir/$project --output_dir $output_dir/$project

for dir in $data_dir/$project/*/; do
    sample=$(basename $dir)
    echo "Processing $sample"
    Rscript scripts/simulate-sce.R \
      --sample_dir $dir \
      --metadata_file $output_dir/$project/single_cell_metadata.tsv \
      --output_dir $output_dir/$project/$sample
done

echo "Creating AnnData files"
Rscript scripts/sce-to-anndata.R --dir $output_dir/$project
python scripts/reformat_anndata.py --dir $output_dir/$project
