#!/bin/bash

# This script runs the neuroblastoma cell type annotation module.
#
# Usage:
#
# ./run-analysis.sh
#
#
# When running in CI or with test data, use:
# testing=1 ./run-analysis.sh
# This will also use the NBAtlas 50K subset to speed up testing.

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define and create directories
data_dir="../../data/current/SCPCP000004"
script_dir="scripts"
ref_dir="references"
scratch_dir="scratch"
results_dir="results"
singler_results_dir="${results_dir}/singler"
mkdir -p $ref_dir
mkdir -p $scratch_dir
mkdir -p $singler_results_dir

# Set up the testing flag
# - If we are testing, we'll use the NBAtlas 50K subset. Otherwise, we'll use the full atlas.
# - We'll also name the NBAtlas reference object files here with the same name as on Mendeley: https://data.mendeley.com/datasets/yhcf6787yp/3
testing=${testing:-0}
if [[ $testing -eq 1 ]]; then
    test_flag="--testing"
    # subset atlas
    nbatlas_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/0a569381-3a0c-4eec-863a-e20544b686ed/file_downloaded"
    nbatlas_seurat="${scratch_dir}/SeuratObj_Share_50kSubset_NBAtlas_v20240130.rds"
else
    test_flag=""
     # full atlas
    nbatlas_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/f5969395-5f6e-4c5d-a61a-5894773d0fee/file_downloaded"
    nbatlas_seurat="${scratch_dir}/seuratObj_NBAtlas_share_v20241203.rds"
fi

nbatlas_sce="${ref_dir}/NBAtlas_sce.rds"
nbatlas_anndata="${ref_dir}/NBAtlas_anndata.h5ad"

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

# define the TumorZoom files
nbatlas_tumor_url="https://data.mendeley.com/public-files/datasets/yhcf6787yp/files/78cad1b4-7425-4073-ba09-362ef73c9ab9/file_downloaded"
nbatlas_tumor_metadata_file="${scratch_dir}/SeuratMeta_Share_TumorZoom_NBAtlas_v20250228.rds"

# Download NBAtlas files
download_file $nbatlas_url $nbatlas_seurat
download_file $nbatlas_tumor_url $nbatlas_tumor_metadata_file

#########################################################
# TODO: This will be going away.
# See: https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/1191
# Determine the ids to retain based on the *full* atlas object
# cell_id_file="${scratch_dir}/nbatlas-cell-ids.txt.gz"
# Rscript ${script_dir}/00a_extract-nbatlas-ids.R \
#     --nbatlas_file "${nbatlas_seurat}" \
#     --cell_id_file "${cell_id_file}"
#########################################################

# Convert the NBAtlas object to SCE
# Temporarily we do not convert to AnnData to save time in CI before we actually get to running scArches
# See: https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/1190
# Rscript ${script_dir}/00b_convert-nbatlas.R \
#    --nbatlas_file "${nbatlas_seurat}" \
#    --tumor_metadata_file "${nbatlas_tumor_metadata_file}" \
#    --cell_id_file "${cell_id_file}" \
#    --sce_file "${nbatlas_sce}" \
#    ${test_flag}
#    # For now, we will not save the AnnData object
#    #--anndata_file "${nbatlas_anndata}"



###################################################################
######################## SingleR annotation #######################
###################################################################


# Train the model with an aggregated NBAtlas reference
singler_model_aggregated="${scratch_dir}/singler-model_nbatlas_aggregated.rds"

# Note can pass in an arbitrary SCE here for sce_file; this is just the first sample in the project
Rscript ${script_dir}/01_train-singler-model.R \
    --nbatlas_sce "${nbatlas_sce}" \
    --sce_file "${data_dir}/SCPCS000101/SCPCL000118_processed.rds" \
    --singler_model_file "${singler_model_aggregated}" \
    --aggregate_reference

# If we are _not_ testing, we'll also train a model with a non-aggregated reference
if [[ $testing == 0 ]]; then
    singler_model_not_aggregated="${scratch_dir}/singler-model_nbatlas_not-aggregated.rds"

    Rscript ${script_dir}/01_train-singler-model.R \
        --nbatlas_sce "${nbatlas_sce}" \
        --sce_file "${data_dir}/SCPCS000101/SCPCL000118_processed.rds" \
        --singler_model_file "${singler_model_not_aggregated}"
fi

# Run SingleR on all samples in the project
for sample_dir in ${data_dir}/SCPCS*; do

    # grab sample id
    sample_id=$(basename $sample_dir)

    # define sample output folder
    sample_results_dir="${singler_results_dir}/${sample_id}"

    for sce_file in $sample_dir/*_processed.rds; do

        # define library ID
        library_id=$(basename $sce_file | sed 's/_processed.rds$//')
        singler_file_basename="${library_id}_singler-annotations.tsv"

        aggr_dir=${sample_results_dir}/aggregated
        mkdir -p $aggr_dir

        # Run with aggregated reference
        singler_output_file="${aggr_dir}/${singler_file_basename}"
        Rscript ${script_dir}/02_classify-singler.R \
            --sce_file "${sce_file}" \
            --singler_model_file "${singler_model_aggregated}" \
            --singler_output_file "${singler_output_file}"

        # If _not_ testing, also run with non-aggregated reference
        if [[ $testing == 0 ]]; then
            nonaggr_dir=${sample_results_dir}/non-aggregated
            mkdir -p $nonaggr_dir
            singler_output_file="${nonaggr_dir}/${singler_file_basename}"

            Rscript ${script_dir}/02_classify-singler.R \
                --sce_file "${sce_file}" \
                --singler_model_file "${singler_model_not_aggregated}" \
                --singler_output_file "{$singler_output_file}"
        fi

    done
done
