#!/bin/bash

# This script runs the module workflow for SCPCP000015 and is meant to be called from ./run-analysis.sh.
# It can be run as follows, optionally specifying the following variables:
#
# threads=<number of threads for inferCNV> seed=<random seed to set for inferCNV> ./run-SCPCP000015.sh

set -euo pipefail

# Ensure script is being run from _one directory up_, which represents the root directory of the module
# This is so `renv` loads properly, without needing to call `renv::load()` in scripts since code run in Docker containers won't use `renv`
module_dir=$(dirname "${BASH_SOURCE[0]}")/..
cd ${module_dir}

# This script processes samples from SCPCP000015
project_id="SCPCP000015"

# Define library to exclude from analysis
exclude_library="SCPCL001111"

# inferCNV settings
threads=${threads:-4}
seed=${seed:-2025}

# Define input directories and files
# These should all be relative to `../`
top_data_dir="../../data/current"
data_dir="../../data/current/${project_id}"
script_dir="scripts"
results_dir="results/${project_id}"
normal_ref_dir="references/normal-references/${project_id}"
merged_sce_file="${top_data_dir}/results/merge-sce/${project_id}/${project_id}_merged.rds"
cell_type_ewings_dir="${top_data_dir}/results/cell-type-ewings/${project_id}"

mkdir -p ${normal_ref_dir}

# Define normal reference files
immune_ref_file="${normal_ref_dir}/ref_immune.rds"
endo_ref_file="${normal_ref_dir}/ref_endo.rds"
endo_immune_ref_file="${normal_ref_dir}/ref_endo-immune.rds"

### Perform the analysis ###

# Build the SCPCP000015 reference files
Rscript ${script_dir}/build-normal-reference/build-reference-SCPCP000015.R \
    --merged_sce_file ${merged_sce_file} \
    --cell_type_ewings_dir ${cell_type_ewings_dir} \
    --reference_immune ${immune_ref_file} \
    --reference_endo ${endo_ref_file} \
    --reference_endo_immune ${endo_immune_ref_file}

# Define all sample ids
sample_ids=$(basename -a ${data_dir}/SCPCS*)

# Run inferCNV on all samples across conditions of interest
for sample_id in $sample_ids; do

    # There is only one library per sample for this project, so we can right away define the SCE
    sce_file=`ls ${data_dir}/${sample_id}/*_processed.rds`

    # skip `exclude_library`
    library_id=$(basename $sce_file | sed 's/_processed.rds$//')
    if [[ $library_id == $exclude_library ]]; then
        continue
    fi

    # Loop over normal references of interest
    # Use reference names (not file names) since these will be used to create directories
    for normal_ref in "ref_immune" "ref_endo" "ref_endo-immune"; do

        # create sample results directory
        sample_results_dir="${results_dir}/${sample_id}/${normal_ref}"
        mkdir -p ${sample_results_dir}

        # define normal reference SCE file
        normal_ref_file="${normal_ref_dir}/${normal_ref}.rds"

        # run inferCNV
        Rscript ${script_dir}/01_run-infercnv.R \
            --sce_file $sce_file \
            --reference_file $normal_ref_file \
            --output_dir $sample_results_dir \
            --threads $threads \
            --seed $seed

        # TODO: Next, we'll run the forthcoming exploratory notebook on each inferCNV result
    done
done
