#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

data_dir="../../data/current"
output_dir="results/simulated"

project=$1

mkdir -p $output_dir/$project

Rscript scripts/permute-metadata.R \
  --metadata_file $data_dir/$project/single_cell_metadata.tsv \
  --output_file $output_dir/$project/single_cell_metadata.tsv

if [ -f $data_dir/$project/bulk_quant.tsv ]; then
  echo "Permuting bulk data"
  Rscript scripts/permute-bulk.R \
    --bulk_file $data_dir/$project/bulk_quant.tsv \
    --output_dir $output_dir/$project
  cp $data_dir/$project/bulk_metadata.tsv $output_dir/$project/bulk_metadata.tsv
else
  echo "Bulk data not found"
fi



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
