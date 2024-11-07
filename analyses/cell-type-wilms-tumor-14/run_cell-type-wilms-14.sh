#!/bin/bash

# This script runs scripts associated with the the cell-type-wilms-tumor-14 module

# set error options
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"
# set CI variable if unset
CI_TESTING=${CI_TESTING:-0}

# set variables for testing
if [ "$CI_TESTING" -eq 0 ]; then
  TEST_FLAG=""
  SCT_FLAG="--run_SCT"
else
  TEST_FLAG="--testing"
  SCT_FLAG=""
fi

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
ref_seurat_sct="${scratch_dir_step}/kidneyatlas_SCT.rdsSeurat"

echo "Downloading reference data"
if [[ ! -e ${ref_h5ad} ]]; then
    ref_url="https://cellgeni.cog.sanger.ac.uk/kidneycellatlas/Fetal_full_v3.h5ad"
    curl -s --show-error --ipv4 -L -o ${ref_h5ad} ${ref_url}
fi

echo "Processing reference data"
Rscript scripts/${step_name}.R \
    --in_fetal_atlas "${ref_h5ad}" \
    --out_fetal_atlas "${ref_seurat}" \
    $SCT_FLAG

echo "Preprocessing data"
## Preprocess data
Rscript scripts/00_preprocessing_rds.R

## Assign anchors

# run specific samples
# Rscript scripts/01_anchor_transfer_seurat.R \
#   --reference "${ref_seurat}" \
#   --metadata "${meta_path}" \
#   --libraries SCPCL000846,SCPCL000847  \
#    $TEST_FLAG

echo "running anchor transfer"
# run all samples
Rscript scripts/01_anchor_transfer_seurat.R \
  --reference "${ref_seurat}" \
  --metadata "${meta_path}" \
  --run_LogNormalize \
  $SCT_FLAG \
  $TEST_FLAG

echo "summarizing results"
Rscript scripts/summary_results.R \
  --metadata "${meta_path}" \
  $TEST_FLAG
