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
consensus_dir="${top_data_dir}/results/cell-type-consensus/${project_id}"

normal_ref_dir="references/normal-references/${project_id}"
mkdir -p ${normal_ref_dir}

# additional input files, again relative to `../`
metadata_file="${data_dir}/single_cell_metadata.tsv"
celltype_reference_group_file="references/reference-cell-groups.tsv"

# Define pooled normal reference files
normal_ref_file="${normal_ref_dir}/all-normal.rds"
immune_ref_file="${normal_ref_dir}/immune.rds"
endo_ref_file="${normal_ref_dir}/endo.rds"

# Define array of references to run across based on $testing
# Use reference names (not file names) since these will be used to create directories below
if [[ $testing -eq 1 ]]; then
    # Use only the all-normal references for testing
    normal_refs=("all-normal")
    test_flag="--testing"
else
    # Use all references for full analysis
    normal_refs=("immune" "endo" "all-normal")
    test_flag=""
fi

# Define all sample ids
sample_ids=$(basename -a ${data_dir}/SCPCS*)

### Perform the analysis ###

# Build the SCPCP000004 normal references
Rscript ${script_dir}/build-normal-reference/build-reference-SCPCP000004.R \
    --merged_sce_file ${merged_sce_file} \
    --celltype_dir ${consensus_dir} \
    --metadata_file ${metadata_file} \
    --celltype_group_file ${celltype_reference_group_file} \
    --reference_normal_rds ${normal_ref_file} \
    --reference_endo_rds ${endo_ref_file} \
    --reference_immune_rds ${immune_ref_file}
