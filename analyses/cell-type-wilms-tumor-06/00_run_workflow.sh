#!/bin/bash

# This script runs the module workflow.
# There are two variables this script will take:
# 1. The TESTING variable controls certain settings to use when this
#  script is being run in CI or more generally being run with test data.
#  Setting TESTING=1 will turn on these settings for test data.
#  This variable is 0 by default.
# 2. The RUN_EXPLORATORY variable controls whether optional exploratory
#  steps that do not directly contribute to final cell type annotations
#  should be run. Setting RUN_EXPLORATORY=1 will run those steps.
#  This variable is 0 by default.
#
# Default usage:
# ./00_run_workflow.sh
#
# Usage in CI:
# TESTING=1 ./00_run_workflow.sh
#
# Usage when running exploratory steps:
# RUN_EXPLORATORY=1 ./00_run_workflow.sh

set -euo pipefail

IS_CI=${TESTING:-0}
RUN_EXPLORATORY=${RUN_EXPLORATORY:-0}
project_id="SCPCP000006"

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define directories
data_dir="../../data/current"
notebook_template_dir="notebook_template"
notebook_output_dir="notebook"

# Download files used for label transfer
# We'll define file names with absolute paths for robustness

# First, download the fetal kidney reference (Stewart et al)
kidney_ref_url="https://datasets.cellxgene.cziscience.com/40ebb8e4-1a25-4a33-b8ff-02d1156e4e9b.rds"
kidney_ref_file="${PWD}/scratch/fetal_kidney.rds"
if [[ ! -f $kidney_ref_file ]]; then
  curl -o $kidney_ref_file $kidney_ref_url
fi
# Second, download the homologs file for gene ID conversion
homologs_url="https://seurat.nygenome.org/azimuth/references/homologs.rds"
homologs_file="${PWD}/scratch/homologs.rds"
if [[ ! -f $homologs_file ]]; then
  curl -o $homologs_file $homologs_url
fi

# Create the fetal references for label transfer (Stewart et al and Cao et al)
Rscript scripts/prepare-fetal-references.R --kidney_ref_file "${kidney_ref_file}"

# Characterize the fetal kidney reference (Stewart et al.)
# This step does not directly contribute to the final annotations
if [[ $RUN_EXPLORATORY -eq 1 ]]; then
  Rscript -e "rmarkdown::render('${notebook_template_dir}/00b_characterize_fetal_kidney_reference_Stewart.Rmd',
      output_format = 'html_document',
      output_file = '00b_characterization_fetal_kidney_reference_Stewart.html',
      output_dir = '${notebook_output_dir}/00-reference',
      params = list(fetal_kidney_path = '${kidney_ref_file}'))"
fi


# Run the label transfer and cluster exploration for all samples in the project
for sample_dir in ${data_dir}/${project_id}/SCPCS*; do
    sample_id=$(basename $sample_dir)

    # define and create sample-specific directories
    # directory for the pre-processed and labeled `Seurat` objects
    results_dir=results/${sample_id}
    # directory for sample-specific notebooks
    sample_notebook_dir=notebook/${sample_id}

    mkdir -p $results_dir
    mkdir -p $sample_notebook_dir

    # Pre-process the data - `Seurat` workflow
    Rscript -e "rmarkdown::render('${notebook_template_dir}/01_seurat-processing.Rmd',
                    params = list(scpca_project_id = '${project_id}', sample_id = '${sample_id}'),
                    output_format = 'html_document',
                    output_file = '01_seurat_processing_${sample_id}.html',
                    output_dir = '${sample_notebook_dir}')"

    # Label transfer from the Cao reference
    Rscript -e "rmarkdown::render('${notebook_template_dir}/02a_label-transfer_fetal_full_reference_Cao.Rmd',
                    params = list(scpca_project_id = '${project_id}', sample_id = '${sample_id}', homologs_file = '${homologs_file}', testing = ${IS_CI}),
                    output_format = 'html_document',
                    output_file = '02a_fetal_all_reference_Cao_${sample_id}.html',
                    output_dir = '${sample_notebook_dir}')"

    # Label transfer from the Stewart reference
    Rscript -e "rmarkdown::render('${notebook_template_dir}/02b_label-transfer_fetal_kidney_reference_Stewart.Rmd',
                    params = list(scpca_project_id = '${project_id}', sample_id = '${sample_id}', homologs_file = '${homologs_file}', testing = ${IS_CI}),
                    output_format = 'html_document',
                    output_file = '02b_fetal_kidney_reference_Stewart_${sample_id}.html',
                    output_dir = '${sample_notebook_dir}')"

    # Cluster exploration
    # This step does not directly contribute to the final annotations
    if [[ $RUN_EXPLORATORY -eq 1 ]]; then
      Rscript -e "rmarkdown::render('${notebook_template_dir}/03_clustering_exploration.Rmd',
                      params = list(scpca_project_id = '${project_id}', sample_id = '${sample_id}', testing = ${IS_CI}),
                      output_format = 'html_document',
                      output_file = '03_clustering_exploration_${sample_id}.html',
                      output_dir = '${sample_notebook_dir}')"
    fi
done


# This step is run here because it must be run both for:
# - scripts/explore-cnv-methods.R (exploratory; calls 06_infercnv.R)
# - 06_infercnv.R (not exploratory)
Rscript scripts/06a_build-geneposition.R

# These steps do not directly contribute to the final annotations
if [[ $RUN_EXPLORATORY -eq 1 ]]; then

  # Run notebook template to explore label transfer and clustering for all samples at once
  for score_threshold in 0.5 0.75 0.85 0.95; do
    Rscript -e "rmarkdown::render('${notebook_output_dir}/04_annotation_Across_Samples_exploration.Rmd',
                    params = list(predicted.score_thr = ${score_threshold}),
                    output_format = 'html_document',
                    output_file = '04_annotation_Across_Samples_exploration_predicted.score_threshold_${score_threshold}.html',
                    output_dir = '${notebook_output_dir}')"
  done

  # Run infercnv and copykat for a selection of samples
  # This script calls scripts/05_copyKAT.R and scripts/06_infercnv.R
  Rscript scripts/explore-cnv-methods.R

fi


# Temporarily, this code is not run in CI
if [[ $IS_CI -eq 0 ]]; then

  # Run infercnv for all samples with HMM i3 and using "both" as the reference
  for sample_dir in ${data_dir}/${project_id}/SCPCS*; do
      sample_id=$(basename $sample_dir)

      # These samples do not have sufficient normal cells to run with a reference in infercnv
      samples_no_reference=("SCPCS000177" "SCPCS000180" "SCPCS000181" "SCPCS000190" "SCPCS000197")

      # Define inferCNV reference set
      if [[ " ${samples_no_reference[*]} " =~ " ${sample_id} " ]]; then
        reference="none"
      else
        reference="both"
      fi

      # don't repeat inference on selection of samples since certain
      #   output files will already exist if exploratory steps were run
      output_file="${results_dir}/${sample_id}/06_infercnv_HMM-i3_${sample_id}_reference-${reference}.rds"
      if [[ ! -f $output_file ]]; then
        Rscript scripts/06_infercnv.R --sample_id $sample_id --reference $reference --HMM i3
      fi
  done

  # Render notebook to make draft annotations
  Rscript -e "rmarkdown::render('${notebook_template_dir}/07_combined_annotation_across_samples_exploration.Rmd',
                          params = list(predicted.celltype.threshold = 0.85, cnv_threshold = 0),
                          output_format = 'html_document',
                          output_file = '07_combined_annotation_across_samples_exploration.html',
                          output_dir = '${notebook_template_dir}')"
fi
