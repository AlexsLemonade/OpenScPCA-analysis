#!/bin/bash

# This script runs the neuroblastoma cell type annotation module.
#
# Usage:
#
# ./run-analysis.sh
#
# When running in CI or with test data, use:
# testing=1 ./run-analysis.sh

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

testing=${testing:-0}

# Define and create directories
script_dir="scripts"
ref_dir="references"
mkdir -p ${ref_dir}

# Define NBAtlas reference object files
nbatlas_seurat="${ref_dir}/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds" # use the file name as recorded in Mendeley
nbatlas_sce="${ref_dir}/NBAtlas_subset_sce.rds"
nbatlas_anndata="${ref_dir}/NBAtlas_subset_anndata.h5ad"
nbatlas_tumor_cells_file="${ref_dir}/NBAtlas_tumor_cells.tsv"

# First, download the NBAtlas Seurat object from: https://data.mendeley.com/datasets/yhcf6787yp/1
# This link downloads the subsetted object with only 50K cells: SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds
nbatlas_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/0a569381-3a0c-4eec-863a-e20544b686ed/file_downloaded"
if [[ ! -f $nbatlas_seurat ]]; thennbatlas_seurat
  curl -o $nbatlas_seurat $nbatlas_url
fi

# Convert the Seurat object to SCE and AnnData formats
Rscript ${script_dir}/00_convert-nbatlas.R
    --seurat_file "${nbatlas_seurat}" \
    --tumor_metadata_file "${nbatlas_tumor_metadata_file}" \
    --sce_file "${nbatlas_sce}" \
    --anndata_file "${nbatlas_anndata}"
