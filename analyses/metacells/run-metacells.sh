#!/bin/bash

set -euo pipefail

# This script is to run scripts in the metacell module, currently in testing

# Run from the script file location
cd "$(dirname "${BASH_SOURCE[0]}")"

data_dir="../../data/current"
results_dir="results"
project="SCPCP000001"
sample="SCPCS000001"
library="SCPCL000001"
test_file="${data_dir}/${project}/${sample}/${library}_processed_rna.h5ad"


# run the SEACells script

python scripts/run-seacells.py \
    $test_file \
    --adata_out ${data_dir}/${project}/${sample}/${library}_seacells.h5ad \
    --model_out ${results_dir}/${project}/${sample}/${library}_seacells_model.pkl \
    --seed 2024
