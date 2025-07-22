#!/bin/bash

# This script runs the neuroblastoma cell type annotation module.
#
# Usage:
#
# ./run-analysis.sh
#
# There are several variables that can be defined when calling this script:
# - testing (Default value: 0)
#   - Use `testing=1` to run with test data; this will also use the NBAtlas 50K subset for efficiency
#   - Example usage: testing=1 ./run-analysis.sh
# - `singler_results_dir` (Default value: "results/singler/")
#   - By default, SingleR results are saved into `results/singler/` in a folder for each sample id.
#     Use this argument to customize the output result directory.
#     The provided value should be _relative_ to the directory this script is stored in.
#   - Example usage: singler_results_dir=relative/other/result/path ./run-analysis.sh
# - `singler_model_file` (Default value: "scratch/singler-model_nbatlas.rds")
#   - By default, the SingleR trained model is saved as `scratch/singler-model_nbatlas.rds`.
#     Use this argument to customize the path/file name which should have an rds file extension.
#     The provided path should be _relative_ to the directory this script is stored in.
#   - Example usage: singler_model_file=path/to/model/file.rds ./run-analysis.sh
# - `aggregate_singler` (Default value: 1)
#   - Use `aggregate_singler=0` to turn off reference aggregation before SingleR model training
#   - Example usage: aggregate_singler=0 ./run-analysis.sh
# - `separate_tumor_singler` (Default value: 0)
#   - Use `separate_tumor_singler=1` to specify that cells in the NBAtlas "tumor zoom" should considered separately from other
#     Neuroendocrine cells during SingleR annotation.
#     By defaut, all Neuroendocrine cells are considered together.
#   - Example usage: separate_tumor_singler=1 ./run-analysis.sh
# - `filter_genes_singler` (Default value: 0)
#   - Use `filter_genes_singler=1` to specify that mitochondrial and ribosomal genes should be removed from NBAtlas
#     before training the SingleR model.
#     By default, these genes are not explicitly removed.
#   - Example usage: filter_genes_singler=1 ./run-analysis.sh
# - `force_convert_nbatlas` (Default value: 0)
#   - This script begins by converting the NBAtlas object to SCE and AnnData formats.
#     By default, if these files exist, the conversion will not be redone.
#     Use this flag to force conversion even if the files already exist.
#   - Example usage: force_convert_nbatlas=1 ./run-analysis.sh
# - `sample_ids` (Default value: "all")
#   - By default, all samples in SCPCP000004 are processed. Specify a set of ids to only run a subset, which
#     may be helpful in certain exploratory runs
#   - Example usage: sample_ids="SCPCS000101 SCPCS000102" ./run-analysis.sh # to run only `SCPCS000101` and `SCPCS000102`
# - `threads` (Default value: 4)
#   - Use `threads=X` to specify a different number of threads
#   - Example usage: threads=2 ./run-analysis.sh # to request 2 threads

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define and create directories
# do this before flag setup since some directory names are used
data_dir="../../data/current/SCPCP000004"
merged_dir="../../data/current/results/merged-sce/SCPCP000004"
script_dir="scripts"
ref_dir="references"
scratch_dir="scratch"
results_dir="results"
mkdir -p $ref_dir
mkdir -p $scratch_dir

###################################################################
######################## Set up variables #########################
###################################################################

# Define argument defaults
testing=${testing:-0} # default is not testing
force_convert_nbatlas=${force_convert_nbatlas:-0} # do not force convert if it already exists
sample_ids=${sample_ids:-"all"} # default is to run all samples
threads=${threads:-4} # default 4 threads

# singler arguments:
singler_results_dir=${singler_results_dir:-"${results_dir}/singler"} # default singler results directory is results/singler
singler_model_file=${singler_model_file:-"${scratch_dir}/singler-model_nbatlas.rds"} # default singler model file is scratch/singler-model_nbatlas.rds
aggregate_singler=${aggregate_singler:-1} # default is to perform aggregation
separate_tumor_singler=${separate_tumor_singler:-0} # default is to _not_ separate tumor cells
filter_genes_singler=${filter_genes_singler:-0} # default is to _not_ filter out genes from NBAtlas


######## Set up singler flags ###########
# Set up singler aggregation
if [[ $aggregate_singler -eq 1 ]]; then
    aggregate_flag="--aggregate_reference"
else
    aggregate_flag=""
fi

# Set up tumor separation setting
if [[ $separate_tumor_singler -eq 1 ]]; then
    separate_tumor_flag="--separate_tumor"
else
    separate_tumor_flag=""
fi

# Set up singler gene filtering
if [[ filter_genes_singler -eq 1 ]]; then
    filter_genes_flag="--filter_genes"
else
    filter_genes_flag=""
fi

####### Set up the testing flag and data ########
# - If we are testing, we'll use the NBAtlas 50K subset. Otherwise, we'll use the full atlas.
# - We'll also name the NBAtlas reference object files here with the same name as on Mendeley:
#   https://data.mendeley.com/datasets/yhcf6787yp/3
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

####### Set up `sample_ids` and check that all sample_ids exist #########
sample_ids=${sample_ids:-"all"}
if [[ $sample_ids == "all" ]]; then
    sample_ids=$(basename -a ${data_dir}/SCPCS*)
else
    for sample_id in ${sample_ids}; do
        if [[ ! -d ${data_dir}/$sample_id ]]; then
            echo "Provided sample_id ${sample_id} isn't available for SCPCP000004."
            echo "Please check for typos or download relevant samples."
            echo "Exiting..."
            exit 1
        fi
    done
fi

###################################################################
######################## Prepare NBAtlas ##########################
###################################################################

echo "Preparing the NBAtlas reference..."

# Define the NBAtlas Seurat and AnnData files
nbatlas_sce="${ref_dir}/NBAtlas_sce.rds"
nbatlas_anndata="${ref_dir}/NBAtlas_anndata.h5ad"
nbatlas_hvg_file="${ref_dir}/NBAtlas_hvgs.txt"

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

# Convert the NBAtlas object to SCE and AnnData formats if needed
if [ ! -f ${nbatlas_sce} ] || [ ! -f ${nbatlas_anndata} ] || [[ ${force_convert_nbatlas} -eq 1 ]]; then
   Rscript ${script_dir}/00_convert-nbatlas.R \
       --nbatlas_file "${nbatlas_seurat}" \
       --tumor_metadata_file "${nbatlas_tumor_metadata_file}" \
       --sce_file "${nbatlas_sce}" \
       --anndata_file "${nbatlas_anndata}" \
       --nbatlas_hvg_file "${nbatlas_hvg_file}" \
       ${test_flag}
fi

###################################################################
######################## SingleR annotation #######################
###################################################################

echo "Training SingleR model..."

# Note can pass in an arbitrary SCE here for sce_file; this is just the first sample in the project
Rscript ${script_dir}/01_train-singler-model.R \
    --nbatlas_sce "${nbatlas_sce}" \
    --sce_file "${data_dir}/SCPCS000101/SCPCL000118_processed.rds" \
    --singler_model_file "${singler_model_file}" \
    --threads $threads \
    ${aggregate_flag} \
    ${separate_tumor_flag} \
    ${filter_genes_flag} \
    ${test_flag}

# Run SingleR on all samples in the project
for sample_id in $sample_ids; do
    echo "Running SingleR on $sample_id..."

    # define sample output folder
    sample_results_dir="${singler_results_dir}/${sample_id}"
    mkdir -p $sample_results_dir

    for sce_file in "${data_dir}/${sample_id}"/*_processed.rds; do

        library_id=$(basename "$sce_file" | sed 's/_processed.rds$//')
        singler_output_tsv="${sample_results_dir}/${library_id}_singler-annotations.tsv"
        singler_output_rds="${sample_results_dir}/${library_id}_singler-results.rds"

        Rscript ${script_dir}/02_classify-singler.R \
            --sce_file "${sce_file}" \
            --singler_model_file "${singler_model_file}" \
            --singler_output_tsv "${singler_output_tsv}" \
            --singler_output_rds "${singler_output_rds}" \
            --threads $threads

    done
done


###################################################################
#################### scANVI/scArches annotation ###################
###################################################################

# Prepare the query merged SCE object for scANVI/scArches

merged_sce_file="${merged_dir}/SCPCP000004_merged.rds"
prepared_anndata_file="${scratch_dir}/SCPCP000004_merged_prepared.h5ad"

Rscript ${script_dir}/03a_prepare-scanvi-query.R \
    --merged_sce_file "${merged_sce_file} \
    --nbatlas_hvg_file "${nbatlas_hvg_file}" \
    --prepared_anndata_file "${prepared_anndata_file}"
