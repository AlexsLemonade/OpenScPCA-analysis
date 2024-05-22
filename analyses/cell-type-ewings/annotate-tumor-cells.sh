#!/bin/bash

####
: '
This workflow will identify tumor cells in a single library.
Tumor cell annotations are obtained by:

- manually identifying cells that express marker genes
- Running CellAssign with a list of tumor marker genes
- Running CopyKAT to identify aneuploid cells

Before running this workflow, you must download both the processed `SingleCellExperiment` objects and `AnnData` objects with `download-data.py`.
Do this for any samples you would like to run through this workflow.

The following arguments are optional and can be used to run this workflow on additional samples (Default sample is "SCPCS000490"):

sample_id: Unique sample ID (name of folder containing libray data)
normal_celltypes: Comma separated list of cell types annotated with either `SingleR` or `CellAssign` to use as a reference list of normal cells.
tumor_celltypes: Comma separated list of cell typs annotated with either `SingleR` or `CellAssign` that are expected to align with tumor cells.
Any cell types used here will be used for comparing to tumor cells annotated in this workflow.

A set of reports summarizing the tumor cell annotations and a TSV file containing annotations from all used methods will be saved to `results/annotate_tumor_cells_output`.

Example of running the workflow with a different sample:

sample_id="SCPCS000491" ./annotate-tumor-cells.sh

'
####

set -euo pipefail

# input variables
sample_id=${sample_id:-"SCPCS000490"}
normal_celltypes=${normal_celltypes:-"Endothelial cells,endothelial cell"}
tumor_celltypes=${tumor_celltypes:-"Pulmonary vascular smooth muscle cells,smooth muscle cell"}

# this script lives in the root of the module directory
# use this to define other default paths
cd $(dirname "$0")
module_directory=$(pwd)

# path to input data folder
data_dir="../../data/current/SCPCP000015"

# define results directories
workflow_results_dir="${module_directory}/results/annotate_tumor_cells_output"
sample_results_dir="${workflow_results_dir}/${sample_id}"

# define output directory for any annotations file
# this is where reference cell tables and annotations files (for InferCNV) will be saved
annotation_dir="$sample_results_dir/annotations"
mkdir -p $annotation_dir

# define paths to notebooks and scripts run in the workflow
notebook_dir="template_notebooks"
scripts_dir="scripts"

# Run preparation scripts to create any references that only need to be created once
# Make gene order file
Rscript $scripts_dir/make-gene-order-file.R

# Run the workflow for each library in the sample directory
for sce in $data_dir/$sample_id/*_processed.rds; do

    # define library ID
    library_id=$(basename $sce | sed 's/_processed.rds$//')

    # Create ref table ---------------------------------------------------

    # define output reference file
    reference_cell_file="$annotation_dir/${library_id}_reference-cells.tsv"

    # Create table with reference cell types
    echo "Starting workflow for $sample_id, $library_id"
    echo "Saving cell references..."
    Rscript $scripts_dir/select-cell-types.R \
      --sce_file "$sce" \
      --normal_cells "${normal_celltypes}" \
      --tumor_cells "${tumor_celltypes}" \
      --output_filename "${reference_cell_file}"

    # Marker gene annotation --------------------------------------------

    # Obtain manual annotations
    echo "Getting manual annotations..."
    Rscript -e "rmarkdown::render('$notebook_dir/01-marker-genes.Rmd', \
          clean = TRUE, \
          output_dir = '$sample_results_dir', \
          output_file = '${library_id}_marker-gene-report.html', \
          params = list(sample_id = '$sample_id', \
                        library_id = '$library_id', \
                        results_dir = '$sample_results_dir', \
                        reference_cell_file = '$reference_cell_file'), \
          envir = new.env()) \
    "


    # InferCNV -------------------------------------------------------

    # define annotations file
    annotations_file="$annotation_dir/${library_id}_normal-annotations.txt"

    # Run InferCNV
    echo "running InferCNV..."
    Rscript $scripts_dir/run-infercnv.Rmd \
      --sce_file "$sce" \
      --annotations_file "$annotations_file" \
      --reference_cell_file "$reference_cell_file" \
      --output_dir "$sample_results_dir/infercnv" \
      --threads 4

    # render infercnv notebook with results
    Rscript -e "rmarkdown::render('$notebook_dir/04-infercnv.Rmd', \
          clean = TRUE, \
          output_dir = '$sample_results_dir', \
          output_file = '${library_id}_infercnv-report.html', \
          params = list(sample_id = '$sample_id', \
                        library_id = '$library_id', \
                        marker_gene_classification_file = '$sample_results_dir/${library_id}_tumor-normal-classifications.tsv', \
                        results_dir = '$sample_results_dir', \
                        infercnv_dir = '$sample_results_dir/infercnv'), \
          envir = new.env()) \
    "

done
