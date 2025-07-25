#!/bin/bash

# This script runs the module workflow for SCPCP000004 and is meant to be called from ./run-analysis.sh.
# It can be run as follows, optionally specifying the following variables:
#
# testing=<0 or 1> threads=<number of threads for inferCNV> seed=<random seed to set for inferCNV> ./run-SCPCP000004.sh

set -euo pipefail

# Ensure script is being run from _one directory up_, which represents the root directory of the module
# This is so `renv` loads properly, without needing to call `renv::load()` in scripts since code run in Docker containers won't use `renv`
module_dir=$(dirname "${BASH_SOURCE[0]}")/..
cd ${module_dir}

# This script processes samples from SCPCP000004
project_id="SCPCP000004"

# Whether to use a subset of references with inferCNV
testing=${testing:-0}

# inferCNV settings
threads=${threads:-4}
seed=${seed:-2025}

# Define input directories and files
# These should all be relative to `../`
top_data_dir="../../data/current"
data_dir="../../data/current/${project_id}"
script_dir="scripts"
notebook_dir="template-notebooks"
results_dir="results/${project_id}"
merged_sce_file="${top_data_dir}/results/merge-sce/${project_id}/${project_id}_merged.rds"

normal_ref_dir="references/normal-references/${project_id}"
mkdir -p ${normal_ref_dir}

# additional input files, again relative to `../`
celltype_reference_group_file="references/reference-cell-groups.tsv"

# Define references to run across based on $testing:
# use only the normal reference for testing, and all references for full analysis
if [[ $testing -eq 1 ]]; then
    normal_refs="normal"
    test_flag="--testing"
else
    normal_refs="immune,endo,normal"
    test_flag=""
fi

# Define groups to include in the `normal` reference
normal_groups="adipocyte,endothelial,epithelial,immune"

# Define low-quality libraries to exclude from references
# context: https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/1234#issuecomment-3113966395
exclude_libraries="SCPCL000124,SCPCL001058"

### Perform the analysis ###

# Build the SCPCP000004 normal references
Rscript ${script_dir}/build-normal-references.R \
    --merged_sce_file ${merged_sce_file} \
    --reference_groups ${normal_refs} \
    --normal_reference_groups ${normal_groups} \
    --reference_celltype_group_file ${celltype_reference_group_file} \
    --reference_output_dir ${normal_ref_dir}
