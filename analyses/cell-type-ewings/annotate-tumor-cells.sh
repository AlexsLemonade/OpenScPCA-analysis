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

# define paths to notebooks and scripts run in the workflow
notebook_dir="template_notebooks"
scripts_dir="scripts"

# define results directories
workflow_results_dir="${module_directory}/results/annotate_tumor_cells_output"
sample_results_dir="${workflow_results_dir}/${sample_id}"

cellassign_results_dir="${sample_results_dir}/cellassign"
mkdir -p $cellassign_results_dir

# define output directory for ref file
ref_dir="${module_directory}/references"
cell_lists_dir="$ref_dir/cell_lists/$sample_id"
mkdir -p $cell_lists_dir

# cellassign refs
all_refs=(
  "${ref_dir}/cellassign_refs/tumor-marker_cellassign.tsv"
  "${ref_dir}/cellassign_refs/visser-all-marker_cellassign.tsv"
  "${ref_dir}/cellassign_refs/panglao-endo-fibro_cellassign.tsv"
)

# Rscript $scripts_dir/generate-cellassign-refs.R

# Run the workflow for each library in the sample directory
for sce in $data_dir/$sample_id/*_processed.rds; do

    # define library ID
    library_id=$(basename $sce | sed 's/_processed.rds$//')

    # path to anndata
    anndata_file="$data_dir/$sample_id/${library_id}_processed_rna.h5ad"

    # define output reference file
    reference_cell_file="$cell_lists_dir/${library_id}_reference-cells.tsv"

    # Create table with reference cell types
    #echo "Starting workflow for $sample_id, $library_id"
    #echo "Saving cell references..."
    #Rscript $scripts_dir/select-cell-types.R \
    #  --sce_file "$sce" \
    #  --normal_cells "${normal_celltypes}" \
    #  --tumor_cells "${tumor_celltypes}" \
    #  --output_filename "${reference_cell_file}"
#
    ## Obtain manual annotations
    #echo "Getting manual annotations..."
    #Rscript -e "rmarkdown::render('$notebook_dir/01-marker-genes.Rmd', \
    #      clean = TRUE, \
    #      output_dir = '$sample_results_dir', \
    #      output_file = '${library_id}_marker-gene-report.html', \
    #      params = list(sample_id = '$sample_id', \
    #                    library_id = '$library_id', \
    #                    results_dir = '$sample_results_dir', \
    #                    reference_cell_file = '$reference_cell_file'), \
    #      envir = new.env()) \
    #"

    # run cell assign for each reference and save predictions file to array
    all_predictions=()

    for i in ${!all_refs[@]}; do

      # grab reference name
      ref_name=$(echo $(basename ${all_refs[$i]}) | sed -e "s/\_cellassign.tsv//")
      echo $ref_name

      # create output predictions file and add to array
      predictions_file="${cellassign_results_dir}/${library_id}_${ref_name}_predictions.tsv"
      all_predictions[$i]=$predictions_file

      # only run cellassign if the predictions file doesn't exist already
      if [ ! -f $predictions_file ]; then
        echo "Running CellAssign for ${library_id} with ${ref_name}"
        python "$scripts_dir/run-cellassign.py" \
          --anndata_file $anndata_file \
          --output_predictions "${predictions_file}" \
          --reference "${all_refs[$i]}" \
          --seed 2024 \
          --threads 4
      fi
    done

    # render report
    Rscript -e "rmarkdown::render('$notebook_dir/02-cellassign.Rmd', \
          clean = TRUE, \
          output_dir = '$sample_results_dir', \
          output_file = '${library_id}_cellassign-report.html', \
          params = list(sample_id = '$sample_id', \
                        library_id = '$library_id', \
                        results_dir = '$sample_results_dir', \
                        tumor_markers_file = '${all_refs[0]}', \ # tumor refs
                        visser_markers_file = '${all_refs[1]}', \ # visser refs
                        marker_gene_classification = '$sample_results_dir/${library_id}_tumor-normal-classifications.tsv', \
                        tumor_marker_predictions = '${predictions_file[0]}', \ # tumor predictions
                        visser_marker_predictions = '${predictions_file[1]}', \ # visser predictions
                        panglao_predictions = '${predictions_file[2]}' \ #panglao predictions
                        ), \
          envir = new.env()) \
    "
done
