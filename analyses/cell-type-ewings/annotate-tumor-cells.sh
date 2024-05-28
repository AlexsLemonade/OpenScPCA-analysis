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
data_dir: Directory containing SCE and AnnData objects to use for a given project.
This is the full path to the project folder.
Default is `../../data/current/SCPCP000015`.
workflow_results_dir: Full path to folder to store outputs from the workflow.
By default, all results will be stored in `results/annotate_tumor_cells_output`.

A set of reports summarizing the tumor cell annotations and a TSV file containing annotations from all used methods will be saved to `results/annotate_tumor_cells_output`.

Example of running the workflow with a different sample:

sample_id="SCPCS000491" ./annotate-tumor-cells.sh

'
####

set -euo pipefail

# this script lives in the root of the module directory
# use this to define other default paths
cd $(dirname "$0")
module_directory=$(pwd)

# input variables
sample_id=${sample_id:-"SCPCS000490"}
normal_celltypes=${normal_celltypes:-"Endothelial cells,endothelial cell"}
tumor_celltypes=${tumor_celltypes:-"Pulmonary vascular smooth muscle cells,smooth muscle cell"}
data_dir=${data_dir:-"../../data/current/SCPCP000015"}
workflow_results_dir=${workflow_results_dir:-"${module_directory}/results/annotate_tumor_cells_output"}
threads=${threads:-4}

echo $workflow_results_dir

# define paths to notebooks and scripts run in the workflow
notebook_dir="template_notebooks"
scripts_dir="scripts"

# define results directories
sample_results_dir="${workflow_results_dir}/${sample_id}"

cellassign_results_dir="${sample_results_dir}/cellassign"
mkdir -p $cellassign_results_dir

# define output directory for any annotations file
# this is where reference cell tables and annotations files (for InferCNV) will be saved
annotation_dir="$sample_results_dir/annotations"
mkdir -p $annotation_dir

# define output directory for refrence matrix files
ref_dir="${module_directory}/references"

# cellassign refs
tumor_only_ref="${ref_dir}/cellassign_refs/tumor-marker_cellassign.tsv"
visser_ref="${ref_dir}/cellassign_refs/visser-all-marker_cellassign.tsv"
panglao_ref="${ref_dir}/cellassign_refs/panglao-endo-fibro_cellassign.tsv"

# Run preparation scripts to create any references that only need to be created once

# generate cell assign refs to use only one time
if [[ ! -f $tumor_only_ref || ! -f $visser_ref || ! -f $panglao_ref ]]; then
  Rscript $scripts_dir/generate-cellassign-refs.R
fi

# Make gene order file if it's not already present
if [ ! -f "$ref_dir/infercnv_refs/Homo_sapiens.GRCh38.104.gene_order.txt" ]; then
  Rscript $scripts_dir/make-gene-order-file.R
fi

# Run the workflow for each library in the sample directory
for sce in $data_dir/$sample_id/*_processed.rds; do

    # define library ID
    library_id=$(basename $sce | sed 's/_processed.rds$//')

    # path to anndata
    anndata_file="$data_dir/$sample_id/${library_id}_processed_rna.h5ad"

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

    # CellAssign ---------------------------------------------------------
    # define output predictions files
    tumor_only_predictions="${cellassign_results_dir}/${library_id}_tumor-marker_predictions.tsv"
    visser_predictions="${cellassign_results_dir}/${library_id}_visser-all-marker_predictions.tsv"
    panglao_predictions="${cellassign_results_dir}/${library_id}_panglao-endo-fibro_predictions.tsv"

    # run cellassign for each reference
    # only run cellassign if the predictions file doesn't exist already
    if [ ! -f $tumor_only_predictions ]; then
      echo "Running CellAssign for ${library_id} with ${tumor_only_ref}"
      python "$scripts_dir/run-cellassign.py" \
        --anndata_file $anndata_file \
        --output_predictions "${tumor_only_predictions}" \
        --reference "${tumor_only_ref}" \
        --seed 2024 \
        --threads $threads
    fi

    if [ ! -f $visser_predictions ]; then
      echo "Running CellAssign for ${library_id} with ${visser_ref}"
      conda run -n openscpca-cell-type-ewings python "$scripts_dir/run-cellassign.py" \
        --anndata_file $anndata_file \
        --output_predictions "${visser_predictions}" \
        --reference "${visser_ref}" \
        --seed 2024 \
        --threads $threads
    fi

    if [ ! -f $panglao_predictions ]; then
      echo "Running CellAssign for ${library_id} with ${panglao_ref}"
      conda run -n openscpca-cell-type-ewings python "$scripts_dir/run-cellassign.py" \
        --anndata_file $anndata_file \
        --output_predictions "${panglao_predictions}" \
        --reference "${panglao_ref}" \
        --seed 2024 \
        --threads $threads
    fi


    # render report
    echo "Rendering CellAssign summary report..."
    Rscript -e "rmarkdown::render('$notebook_dir/02-cellassign.Rmd', \
          clean = TRUE, \
          output_dir = '$sample_results_dir', \
          output_file = '${library_id}_cellassign-report.html', \
          params = list(sample_id = '$sample_id', \
                        library_id = '$library_id', \
                        results_dir = '$sample_results_dir', \
                        marker_gene_classification = '$sample_results_dir/${library_id}_tumor-normal-classifications.tsv', \
                        tumor_marker_predictions = '$tumor_only_predictions', \
                        visser_marker_predictions = '$visser_predictions', \
                        panglao_predictions = '$panglao_predictions'), \
          envir = new.env()) \
    "

    # CopyKAT ----------------------------------------------------------
    if [ ! -f "$sample_results_dir/copykat/no_reference/${library_id}_final-copykat.rds" ]; then
      echo "Running CopyKAT with no reference..."
      Rscript $scripts_dir/run-copykat.R \
        --sce_file "$sce" \
        --results_dir "$sample_results_dir/copykat/no_reference" \
        --threads $threads
    fi

    if [ ! -f "$sample_results_dir/copykat/with_reference/${library_id}_final-copykat.rds" ]; then
      echo "Running CopyKAT with a reference..."
      Rscript $scripts_dir/run-copykat.R \
        --sce_file "$sce" \
        --reference_cell_file "$reference_cell_file" \
        --results_dir "$sample_results_dir/copykat/with_reference" \
        --threads $threads
    fi

    echo "Rendering CopyKAT report..."
    Rscript -e "rmarkdown::render('$notebook_dir/03-copykat.Rmd', \
          clean = TRUE, \
          output_dir = '$sample_results_dir', \
          output_file = '${library_id}_copykat-report.html', \
          params = list(sample_id = '$sample_id', \
                        library_id = '$library_id', \
                        marker_gene_classification = '$sample_results_dir/${library_id}_tumor-normal-classifications.tsv', \
                        reference_cell_file = '$reference_cell_file', \
                        no_ref_copykat_results = '$sample_results_dir/copykat/no_reference', \
                        with_ref_copykat_results = '$sample_results_dir/copykat/with_reference', \
                        results_dir = '$sample_results_dir'), \
          envir = new.env()) \
    "


    # InferCNV -------------------------------------------------------

    # define annotations file
    annotations_file="$annotation_dir/${library_id}_normal-annotations.txt"

    # Run InferCNV
    if [ ! -f "$sample_results_dir/infercnv/${library_id}_cnv-obj.rds" ]; then
      echo "running InferCNV..."
      Rscript $scripts_dir/run-infercnv.R \
        --sce_file "$sce" \
        --annotations_file "$annotations_file" \
        --reference_cell_file "$reference_cell_file" \
        --output_dir "$sample_results_dir/infercnv" \
        --threads 4
    fi

    # render infercnv notebook with results
    Rscript -e "rmarkdown::render('$notebook_dir/04-infercnv.Rmd', \
          clean = TRUE, \
          output_dir = '$sample_results_dir', \
          output_file = '${library_id}_infercnv-report.html', \
          params = list(sample_id = '$sample_id', \
                        library_id = '$library_id', \
                        marker_gene_classification = '$sample_results_dir/${library_id}_tumor-normal-classifications.tsv', \
                        results_dir = '$sample_results_dir', \
                        infercnv_dir = '$sample_results_dir/infercnv'), \
          envir = new.env()) \
    "

done
