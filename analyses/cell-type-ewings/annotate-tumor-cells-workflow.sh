#!/usr/bin/env bash

####
: '
This workflow will identify tumor cells in a single library.
Tumor cell annotations are obtained by:

- manually identifying cells that express marker genes
- Running CellAssign with a list of tumor marker genes
- Running CopyKAT to identify aneuploid cells

Before running this workflow, you must download both the processed `SingleCellExperiment` objects and `AnnData` objects with `download-data.py`.
Do this for any samples you would like to run through this workflow.

Then create a tsv file with the following columns for each sample that you would like to run through the worfklow:

sample_id: Unique sample ID (name of folder containing libray data)
library_id: Unique library ID (prefix of data files)
normal_celltypes: Comma separated list of cell types annotated with either `SingleR` or `CellAssign` to use as a reference list of normal cells.
tumor_celltypes: Comma separated list of cell typs annotated with either `SingleR` or `CellAssign` that are expected to align with tumor cells.
Any cell types used here will be used for comparing to tumor cells annotated in this workflow.
If you do not wish to provide any tumor cell types, leave blank.

The path to this file should be provided when running the workflow using the --sample_metadata argument.

Optionally, you can include the `data_dir` argument. Use this if you want to use a different release other than `current` after downloading data.

A set of reports summarizing the tumor cell annotations and a TSV file containing annotations from all used methods will be saved to the `results_dir`.
The default `results_dir` is `results/annotate_tumor_cells_output`.

Example of running the workflow:

./annotate-tumor-cells-workflow.sh \
  --sample_metadata "sample_metadata.tsv"

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

# results and refs
ref_dir="${module_directory}/references"
workflow_results_dir="${module_directory}/results/annotate_tumor_cells_output"

# define paths to notebooks and scripts run in the workflow
notebook_dir="template_notebooks"
scripts_dir="scripts"

# Run the workflow for each library in the sample directory
for sce in $data_dir/$sample_id/*_processed.rds; do

    # define library ID
    library_id=$(basename $sce | sed 's/_processed.rds$//')

    # define output directory for ref file
    cell_lists_dir="$ref_dir/cell_lists/$sample_id"
    # define output reference file
    reference_cell_file="$cell_lists_dir/${library_id}_reference-cells.tsv"

    # directory to save all results for sample
    sample_results_dir="${workflow_results_dir}/${sample_id}"

    # Create table with reference cell types
    echo "Starting workflow for $sample_id, $library_id"
    echo "Saving cell references..."
    Rscript $scripts_dir/select-cell-types.R \
      --sce_file "$sce" \
      --normal_cells "${normal_celltypes}" \
      --tumor_cells "${tumor_celltypes}" \
      --output_filename "${reference_cell_file}"

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
done
