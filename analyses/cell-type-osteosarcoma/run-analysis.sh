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

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

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
scpca_sce_file=$(find "${data_dir}/${project_id}" -name "*_processed.rds" | sort | head -n 1)

# Define argument defaults
testing=${testing:-0} # default is not testing
threads=${threads:-4} # default 4 threads

# Define argument defaults
singler_model_file=${singler_model_file:-"${scratch_dir}/singler-model_osteocar.rds"} # default singler model file is scratch/singler-model_osteocar.rds


# Download OsteoCar reference ----------------------------------------------
figshare_mets_url="https://api.figshare.com/v2/file/download/65376111"
figshare_prim_url="https://api.figshare.com/v2/file/download/65375988"
figshare_xeno_mets_url="https://api.figshare.com/v2/file/download/65375973"
figshare_xeno_prim_url="https://api.figshare.com/v2/file/download/65375976"
ref_urls=(
  $figshare_mets_url
  $figshare_prim_url
  $figshare_xeno_mets_url
  $figshare_xeno_prim_url
)
ref_prefixes=(
  "${ref_dir}/patient_mets"
  "${ref_dir}/patient_prim"
  "${ref_dir}/xeno_mets"
  "${ref_dir}/xeno_prim"
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

# Create harmonized, aggregated SingleR references -------------------------
# The primary and metastasis references are combined into one reference with
# harmonized labels. This is currently only done for the patient references.
Rscript "${script_dir}/02-create-reference.R" \
  --prim_sce_file "${ref_dir}/patient_prim_sce.rds" \
  --mets_sce_file "${ref_dir}/patient_mets_sce.rds" \
  --scpca_sce_file "${scpca_sce_file}" \
  --output_ref_file "${ref_dir}/patient_harmonized_aggr_ref.rds" \
  --threads $threads

# Train SingleR model ------------------------------------------------------
# The ScPCA SCE defines the genes the model expects in query data
Rscript "${script_dir}/03-train-singler-model.R" \
  --ref_file "${ref_dir}/patient_harmonized_aggr_ref.rds" \
  --sce_file "${scpca_sce_file}" \
  --singler_model_file "${singler_model_file}" \
  --threads $threads
