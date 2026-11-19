#!/bin/bash
 
# This script runs the osteosarcoma cell type annotation module.
#
# Usage:
#
# ./run-analysis.sh
#
# There are several variables that can be defined when calling this script:
# - testing (Default value: 0)
#   - Use `testing=1` to run with test data
#   - Example usage: testing=1 ./run-analysis.sh
# - `threads` (Default value: 4)
#   - Use `threads=X` to specify a different number of threads
#   - Example usage: threads=2 ./run-analysis.sh # to request 2 threads
# - `sample_type` (Default value: patient)
#   - Which OsteoCAR references to use, either `patient` or `xeno`
#   - Only the references for this sample type are downloaded and processed
#   - Example usage: sample_type=xeno ./run-analysis.sh
 
set -euo pipefail
 
# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}
 
# Define argument defaults
testing=${testing:-0}               # default is not testing
threads=${threads:-4}               # default 4 threads
sample_type=${sample_type:-patient} # default is the patient references
 
# Check sample_type
if [[ $sample_type != "patient" && $sample_type != "xeno" ]]; then
  echo "sample_type must be either 'patient' or 'xeno'"
  exit 1
fi
 
# Define and create directories
script_dir="scripts"
ref_dir="references"
scratch_dir="scratch"
mkdir -p $ref_dir
mkdir -p $scratch_dir
 
# ScPCA data for the osteosarcoma project
# Any processed SCE works here: all ScPCA libraries share the same genes,
# so this one is only used to define the gene symbols the query data will have
data_dir="../../data/current"
project_id="SCPCP000017"
scpca_sce_file="${data_dir}/${project_id}/SCPCS000246/SCPCL000338_processed.rds"
 
# Output files, named by sample type
ref_prefix="${ref_dir}/${sample_type}"
harmonized_ref_file="${ref_prefix}_harmonized_aggr_ref.rds"
singler_model_file=${singler_model_file:-"${scratch_dir}/${sample_type}_singler-model_osteocar.rds"}
 
 
# Download OsteoCar reference ----------------------------------------------
# Only the primary and metastasis references for the chosen sample type are used
if [[ $sample_type == "patient" ]]; then
  figshare_prim_url="https://api.figshare.com/v2/file/download/65375988"
  figshare_mets_url="https://api.figshare.com/v2/file/download/65376111"
else
  figshare_prim_url="https://api.figshare.com/v2/file/download/65375976"
  figshare_mets_url="https://api.figshare.com/v2/file/download/65375973"
fi

ref_urls=(
  $figshare_prim_url
  $figshare_mets_url
)
ref_prefixes=(
  "${ref_prefix}_prim"
  "${ref_prefix}_mets"
)
 
# First, download the object with a helper function
# This function takes two arguments in order, the URL and the filename to save to
download_file() {
  local url="$1"
  local file_name="$2"
  if [[ ! -f $file_name ]]; then
    curl -Lo $file_name $url
  fi
}
 
# Download OsteoCar reference files
for i in "${!ref_urls[@]}"; do
  download_file "${ref_urls[$i]}" "${ref_prefixes[$i]}.qs2"
done
 
# Convert to SCE and AnnData -----------------------------------------------
for i in "${!ref_prefixes[@]}"; do
  Rscript "${script_dir}/01-convert-osteocar.R" \
    --input_ref_file "${ref_prefixes[$i]}.qs2" \
    --output_sce_file "${ref_prefixes[$i]}_sce.rds" \
    --output_anndata_file "${ref_prefixes[$i]}_adata.h5ad"
done
 
# Create harmonized, aggregated SingleR reference --------------------------
# The primary and metastasis references are combined into one reference with
# harmonized labels.
Rscript "${script_dir}/02-create-reference.R" \
  --prim_sce_file "${ref_prefix}_prim_sce.rds" \
  --mets_sce_file "${ref_prefix}_mets_sce.rds" \
  --output_ref_file "${harmonized_ref_file}" \
  --threads $threads
 
# Train SingleR model ------------------------------------------------------
# The ScPCA SCE defines the genes the model expects in query data
Rscript "${script_dir}/03-train-singler-model.R" \
  --ref_file "${harmonized_ref_file}" \
  --sce_file "${scpca_sce_file}" \
  --singler_model_file "${singler_model_file}" \
  --threads $threads
