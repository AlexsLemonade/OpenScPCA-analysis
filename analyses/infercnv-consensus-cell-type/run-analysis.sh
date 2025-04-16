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

# By default, use 4 cores
threads=4

# Define directories and file paths
data_dir="../../data/current"
script_dir="scripts"
results_dir="results"
ref_dir="references"
infercnv_ref_dir="${ref_dir}/normal-references"

mkdir -p ${infercnv_ref_dir}

# Create the gene order file for input to inferCNV
Rscript ${script_dir}/00-make-gene-order-file.R

##### Analysis for SCPCP000015 #####

project_id="SCPCP000015"
project_data_dir="${data_dir}/${project_id}/"
project_results_dir="${results_dir}/${project_id}/"
mkdir -p ${project_results_dir}

# Define library to exclude from inferCNV
exclude_library="SCPCL001111"

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

# Run inferCNV across libraries
sample_ids=$(basename -a ${project_data_dir}/SCPCS*)
for sample_id in $sample_ids; do

    # top-level sample results directory
    sample_results_dir="${project_results_dir}/${sample_id}"

    for sce_file in ${project_data_dir}/${sample_id}/*_processed.rds; do
        # skip ${exclude_library}
        library_id=$(basename $sce_file | sed 's/_processed.rds$//')
        if [[ $library_id == $exclude_library ]]; then
            continue
        fi

        # run infercnv using both references
        Rscript ${script_dir}/01_run-infercnv.R \
            --sce_file $sce_file \
            --reference_file $immune_ref_file \
            --output_dir $sample_results_dir/ref-all-immune \
            --hmm_model "i6" \
            --threads $threads

        Rscript ${script_dir}/01_run-infercnv.R \
            --sce_file $sce_file \
            --reference_file $immune_subset_ref_file \
            --output_dir $sample_results_dir/ref-subset-immune \
            --hmm_model "i6" \
            --threads $threads

    done
done

