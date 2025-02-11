#!/bin/bash

# This script runs the module workflow.
# There are two variables this script will take:
# 1. The TESTING variable controls certain settings to use when this
#  script is being run in CI or more generally being run with test data.
#  Setting TESTING=1 will turn on these settings for test data.
#  This variable is 0 by default.
# 2. The RUN_EXPLORATORY variable controls whether optional exploratory
#  steps that generate notebooks stored in the repository but do not
#  directly contribute to final cell type annotations should be run.
#  Setting RUN_EXPLORATORY=1 will run those steps.
#  This variable is 0 by default.
# 3. The THREADS variable controls how many cores are used for inference
#  with copyKAT, which is an exploratory step in the workflow controlled by RUN_CNV_EXPLORATORY=1.
#  The variable is 32 by default.
#
# Default usage:
# ./00_run_workflow.sh
#
# Usage in CI:
# TESTING=1 ./00_run_workflow.sh
#
# Usage when running exploratory steps:
# This mode should be used to regenerate all notebooks in the repository.
# RUN_EXPLORATORY=1 ./00_run_workflow.sh

set -euo pipefail

TESTING=${TESTING:-0}
RUN_EXPLORATORY=${RUN_EXPLORATORY:-0}
THREADS=${THREADS:-32}
project_id="SCPCP000006"

# Define the `predicted_celltype_threshold` which is used to identify cells with a high confident label transfer from the fetal kidney reference.
# This score is used in the script `06_infercnv.R` and in the notebook `07_combined_annotation_across_samples.Rmd`
predicted_celltype_threshold=0.75

# Define the `cnv_threshold_low/high` which are used to identify cells with no CNV (`cnv_score` <= `cnv_threshold_low`) or with CNV (`cnv_score` > `cnv_threshold_high`)
cnv_threshold_low=0
cnv_threshold_high=2

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define directories
data_dir="../../data/current"
notebook_template_dir="notebook_template"
notebook_dir="notebook"

# Define test data string to use with 06_infercnv.R
if [[ $TESTING -eq 1 ]]; then
  test_string="--testing"
else
  test_string=""
fi

# Download files used for label transfer
# We'll define file names with absolute paths for robustness


# Prepare gene file for inferCNV
Rscript scripts/06a_build-geneposition.R

# Run infercnv for all samples with HMM i3 and using "both" as the reference, where possible
for sample_dir in ${data_dir}/${project_id}/SCPCS*; do
    sample_id=$(basename $sample_dir)
    results_dir=results/${sample_id}

    # These samples do not have sufficient normal cells to run with a reference in infercnv
    samples_no_reference=("SCPCS000177" "SCPCS000180" "SCPCS000181" "SCPCS000190" "SCPCS000197")

    # Define inferCNV reference set
    if [[ " ${samples_no_reference[*]} " =~ " ${sample_id} " ]]; then
      reference="none"
    else
      reference="both"
    fi

    # Run inferCNV
    Rscript scripts/06_infercnv.R --sample_id $sample_id --reference $reference --HMM i3 ${test_string} --score_threshold ${predicted_celltype_threshold}
done


# Render notebook to make draft annotations
Rscript -e "rmarkdown::render('${notebook_dir}/07_combined_annotation_across_samples_exploration.Rmd',
                        params = list(predicted.celltype.threshold = ${predicted_celltype_threshold}, cnv_threshold_low = ${cnv_threshold_low}, cnv_threshold_high = ${cnv_threshold_high}, testing = ${TESTING}),
                        output_format = 'html_document',
                        output_file = '07_combined_annotation_across_samples_exploration.html',
                        output_dir = '${notebook_dir}')"
