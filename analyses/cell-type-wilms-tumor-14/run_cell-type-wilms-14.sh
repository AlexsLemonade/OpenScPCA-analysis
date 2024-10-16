#!/bin/bash

# This script runs scripts associated with the the cell-type-wilms-tumor-14 module 

# set error options
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
# set CI variable if unset
CI_TESTING=${CI_TESTING:-0}

scratch_dir="scratch"
results_dir="results"
data_dir="../../data/current/SCPCP000014"

meta_path="${data_dir}/single_cell_metadata.tsv"



## Preprocess reference single-cell datasets (fetal atlas)
step_name="00_preprocess_reference"
scratch_dir_step="${scratch_dir}/${step_name}" && mkdir -p ${scratch_dir_step}

# Download and process reference data
ref_h5ad="${scratch_dir_step}/Fetal_full_v3.h5ad" 
ref_seurat="${scratch_dir_step}/kidneyatlas.rdsSeurat" 

if [[ ! -e ${ref_h5ad} ]]; then
    ref_url="https://cellgeni.cog.sanger.ac.uk/kidneycellatlas/Fetal_full_v3.h5ad"
    curl -o ${ref_h5ad} ${ref_url}
fi

Rscript scripts/${step_name}.R \
    --in_fetal_atlas "${ref_h5ad}" \
    --out_fetal_atlas "${ref_seurat}"

## Preprocess data
Rscript scripts/00_preprocessing_rds.R

## Assign anchors 
if [ "$CI_TESTING" -eq 0 ]; then
  TEST_FLAG=""
else
  TEST_FLAG="--testing"
fi


# run specific samples
# Rscript scripts/01_anchor_transfer_seurat.R \
#   --reference "${ref_seurat}" \
#   --metadata "${meta_path}" \
#   --libraries SCPCL000846,SCPCL000847  \
#    $TEST_FLAG

# run all samples
Rscript scripts/01_anchor_transfer_seurat.R \
  --reference "${ref_seurat}" \
  --metadata "${meta_path}" \
  $TEST_FLAG