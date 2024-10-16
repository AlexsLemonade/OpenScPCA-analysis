#!/bin/bash

# This script runs the module workflow.
#
# USAGE:
# bash 00_run_workflow.sh
#
# USAGE in CI:
# TESTING=1 bash 00_run_workflow.sh

set -euo pipefail

IS_CI=${TESTING:-0}
project_id="SCPCP000006"

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define directories
data_dir="../../data/current"
notebook_template_dir="notebook_template"
notebook_output_dir="notebook"


# Download and create the fetal kidney reference (Stewart et al)
Rscript scripts/download-and-create-fetal-kidney-ref.R

# Build the gene position file reference for infercnv
Rscript scripts/06a_build-geneposition.R

# Characterize the fetal kidney reference (Stewart et al.)
Rscript -e "rmarkdown::render('${notebook_template_dir}/00b_characterize_fetal_kidney_reference_Stewart.Rmd',
    output_format = 'html_document',
    output_file = '00b_characterization_fetal_kidney_reference_Stewart.html',
    output_dir = '${notebook_output_dir}/00-reference')"


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

    # Temporarily this code is not run in CI.
    if [[ $IS_CI -eq 0 ]]; then

        # Label transfer from the Cao reference using Azimuth
        Rscript -e "rmarkdown::render('${notebook_template_dir}/02a_label-transfer_fetal_full_reference_Cao.Rmd',
                        params = list(scpca_project_id = '${project_id}', sample_id = '${sample_id}'),
                        output_format = 'html_document',
                        output_file = '02a_fetal_all_reference_Cao_${sample_id}.html',
                        output_dir = '${sample_notebook_dir}')"

        # Label transfer from the Stewart reference using Seurat
        Rscript -e "rmarkdown::render('${notebook_template_dir}/02a_label-transfer_fetal_full_reference_Stewart.Rmd',
                        params = list(scpca_project_id = '${project_id}', sample_id = '${sample_id}'),
                        output_format = 'html_document',
                        output_file = '02a_fetal_all_reference_Stewart_${sample_id}.html',
                        output_dir = '${sample_notebook_dir}')"

        # Cluster exploration
        Rscript -e "rmarkdown::render('${notebook_template_dir}/03_clustering_exploration.Rmd',
                        params = list(scpca_project_id = '${project_id}', sample_id = '${sample_id}'),
                        output_format = 'html_document',
                        output_file = '03_clustering_exploration_${sample_id}.html',
                        output_dir = '${sample_notebook_dir}')"
    fi
done

# Temporarily this code is not run in CI.
if [[ $IS_CI -eq 0 ]]; then

  # Run notebook template to explore label transfer and clustering for all samples at once
  Rscript -e "rmarkdown::render('${notebook_output_dir}/04_annotation_Across_Samples_exploration.Rmd',
                  output_format = 'html_document',
                  output_file = '04_annotation_Across_Samples_exploration.html',
                  output_dir = ${notebook_output_dir})"

  # Run infercnv and copykat for a selection of samples
  # This script calls scripts/05_copyKAT.R and scripts/06_infercnv.R
  Rscript scripts/explore-cnv-methods.R

fi
