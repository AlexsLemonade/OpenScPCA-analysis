#!/bin/bash

# This script runs the module workflow.
#
# Usage:
#
# ./run-analysis.sh

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define directories and file paths
data_dir="../../data/current"
script_dir="scripts"
results_dir="results"
infercnv_ref_dir="${results_dir}/normal-references"

mkdir -p ${infercnv_ref_dir}

##### Analysis for SCPCP000015 #####

# Define input files for scripts
merged_sce_file="${data_dir}/results/merge-sce/SCPCP000015/SCPCP000015_merged.rds"
cell_type_ewings_dir="${data_dir}/results/cell-type-ewings/SCPCP000015"

# Define normal reference files
ewings_ref_dir="${infercnv_ref_dir}/SCPCP000015"
mkdir -p ${ewings_ref_dir}
immune_ref_file="${ewings_ref_dir}/ref-all-immune.rds"
immune_subset_ref_file="${ewings_ref_dir}/ref-subset-immune.rds"

# Build the SCPCP000015 reference files
Rscript ${script_dir}/build-normal-reference/build-reference-SCPCP000015.R \
    --merged_sce_file ${merged_sce_file} \
    --cell_type_ewings_dir ${cell_type_ewings_dir} \
    --reference_immune ${immune_ref_file} \
    --reference_immune_subset ${immune_subset_ref_file}
