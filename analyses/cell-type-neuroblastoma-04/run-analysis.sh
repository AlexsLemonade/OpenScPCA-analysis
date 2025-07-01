#!/bin/bash

# This script runs the neuroblastoma cell type annotation module.
#
# Usage:
#
# ./run-analysis.sh
#
# By default, this will run with the full NBAtlas dataset. To use the 50k subset instead, use:
# nbatlas_version="subset" ./run-analysis.sh
#
# When running in CI or with test data, use:
# testing=1 ./run-analysis.sh

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define the nbatlas_version
nbatlas_version=${nbatlas_version:-"full"}

# Set up the testing flag
testing=${testing:-0}
if [[ $testing -eq 1 ]]; then
    test_flag="--testing"
else
    test_flag=""
fi

# Define and create directories
data_dir="../../data/current/SCPCP000004"
script_dir="scripts"
ref_dir="references"
scratch_dir="scratch"
mkdir -p $ref_dir
mkdir -p $scratch_dir

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

nbatlas_sce="${ref_dir}/NBAtlas_${nbatlas_version}_sce.rds"
nbatlas_anndata="${ref_dir}/NBAtlas_${nbatlas_version}_anndata.h5ad"

###################################################################
######################## Prepare NBAtlas ##########################
###################################################################

# First, download the NBAtlas Seurat objects from Mendeley with a helper function
# This function takes two arguments in order, the URL and the filename to save to
download_file() {
  local url="$1"
  local file_name="$2"
  if [[ ! -f $file_name ]]; then
    curl -L -o $file_name $url
  fi
}

# Define URLs; the first 2 will always be downloaded
nbatlas_full_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/f5969395-5f6e-4c5d-a61a-5894773d0fee/file_downloaded"
nbatlas_tumor_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/78cad1b4-7425-4073-ba09-362ef73c9ab9/file_downloaded"
nbatlas_subset_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/0a569381-3a0c-4eec-863a-e20544b686ed/file_downloaded"

# Define NBAtlas reference object files with the same name as on Mendeley: https://data.mendeley.com/datasets/yhcf6787yp/3
nbatlas_full_seurat="${scratch_dir}/seuratObj_NBAtlas_share_v20241203.rds" # this is needed regardless of which nbatlas_version was specified since we need to subset ids
nbatlas_tumor_metadata_file="${scratch_dir}/SeuratMeta_Share_TumorZoom_NBAtlas_v20250228.rds"

# Define which NBAtlas version we're actually converting, and download subset if it's that one
if [[ $nbatlas_version == "full" ]]; then
  nbatlas_seurat=${nbatlas_full_seurat}
else
  nbatlas_seurat="${scratch_dir}/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds"
  download_file $nbatlas_subset_url $nbatlas_seurat
fi

# grab remaining downloads
download_file $nbatlas_full_url $nbatlas_full_seurat
download_file $nbatlas_tumor_url $nbatlas_tumor_metadata_file


# Determine the ids to retain based on the *full* atlas object
cell_id_file="${scratch_dir}/nbatlas-cell-ids.txt.gz"
Rscript ${script_dir}/00a_extract-nbatlas-ids.R \
    --nbatlas_file "${nbatlas_full_seurat}" \
    --cell_id_file "${cell_id_file}"

# Convert the NBAtlas object to SCE (and in the future, also to AnnData)
Rscript ${script_dir}/00b_convert-nbatlas.R \
   --nbatlas_file "${nbatlas_seurat}" \
   --tumor_metadata_file "${nbatlas_tumor_metadata_file}" \
   --cell_id_file "${cell_id_file}" \
   --sce_file "${nbatlas_sce}" \
   ${test_flag}
   # For now, we will not save the AnnData object
   #--anndata_file "${nbatlas_anndata}"



###################################################################
######################## SingleR annotation #######################
###################################################################

echo "Training SingleR model"

# Train the SingleR model with an aggregated reference
singler_model_file="${scratch_dir}/singler-model_aggregated_nbatlas-${nbatlas_version}.rds"

# can pass in an arbitrary SCE here for sce_file; this is just the first sample in the project
Rscript ${script_dir}/01_train-singler-model.R \
    --nbatlas_sce "${nbatlas_sce}" \
    --sce_file "${data_dir}/SCPCS000101/SCPCL000118_processed.rds" \
    --singler_model_file "${singler_model_file}" \
    --aggregate_reference
