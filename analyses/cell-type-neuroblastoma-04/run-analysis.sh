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

# Set up the testing flag
testing=${testing:-0}
if [[ $testing -eq 1 ]]; then
    test_flag="--testing"
else
    test_flag=""
fi

# Define and create directories
script_dir="scripts"
ref_dir="references"
scratch_dir="scratch"
mkdir -p $ref_dir
mkdir -p $scratch_dir

# Define NBAtlas reference object files
# Ensure local files have the same name as on Mendeley: https://data.mendeley.com/datasets/yhcf6787yp/3
nbatlas_subset_seurat="${scratch_dir}/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds"
nbatlas_full_seurat="${scratch_dir}/seuratObj_NBAtlas_share_v20241203.rds"
nbatlas_tumor_metadata_file="${scratch_dir}/SeuratMeta_Share_TumorZoom_NBAtlas_v20250228.rds"

nbatlas_subset_sce="${ref_dir}/NBAtlas_subset_sce.rds"
nbatlas_subset_anndata="${ref_dir}/NBAtlas_subset_anndata.h5ad"
nbatlas_full_sce="${ref_dir}/NBAtlas_full_sce.rds"
nbatlas_full_anndata="${ref_dir}/NBAtlas_full_anndata.h5ad"

# First, download the NBAtlas Seurat objects from Mendeley with a helper function
# This function takes two arguments in order, the URL and the filename to save to
download_file() {
  local url="$1"
  local file_name="$2"
  if [[ ! -f $file_name ]]; then
    curl -L -o $file_name $url
  fi
}

nbatlas_subset_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/0a569381-3a0c-4eec-863a-e20544b686ed/file_downloaded"
nbatlas_full_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/f5969395-5f6e-4c5d-a61a-5894773d0fee/file_downloaded"
nbatlas_tumor_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/78cad1b4-7425-4073-ba09-362ef73c9ab9/file_downloaded"

download_file $nbatlas_subset_url $nbatlas_subset_seurat
download_file $nbatlas_full_url $nbatlas_full_seurat
download_file $nbatlas_tumor_url $nbatlas_tumor_metadata_file

# Convert the NBAtlas object to SCE and AnnData formats, for both the full and subsetted versions
Rscript ${script_dir}/00_convert-nbatlas.R \
    --seurat_file "${nbatlas_subset_seurat}" \
    --tumor_metadata_file "${nbatlas_tumor_metadata_file}" \
    --sce_file "${nbatlas_subset_sce}" \
    --anndata_file "${nbatlas_subset_anndata}" \
    ${test_flag}

Rscript ${script_dir}/00_convert-nbatlas.R \
    --seurat_file "${nbatlas_full_seurat}" \
    --tumor_metadata_file "${nbatlas_tumor_metadata_file}" \
    --sce_file "${nbatlas_full_sce}" \
    --anndata_file "${nbatlas_full_anndata}" \
    ${test_flag}
