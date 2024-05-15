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

To run the workflow on an individual sample provide the following arguments:

sample_id: Unique sample ID (name of folder containing libray data)
library_id: Unique library ID (prefix of data files)
normal_celltypes: Comma separated list of cell types annotated with either `SingleR` or `CellAssign` to use as a reference list of normal cells.
tumor_celltypes: Comma separated list of cell typs annotated with either `SingleR` or `CellAssign` that are expected to align with tumor cells.
Any cell types used here will be used for comparing to tumor cells annotated in this workflow.

Optionally, you can include the `data_dir` column. Use this if you want to use a different release other than `current` after downloading data.

A set of reports summarizing the tumor cell annotations and a TSV file containing annotations from all used methods will be saved to the `results_dir`.
The default `results_dir` is `results/annotate_tumor_cells_output/<sample_id>`.

Example of running the workflow:

./annotate-tumor-cells-workflow.sh \
  --sample_id "SCPCS000490" \
  --library_id "SCPCL000822" \
  --normal_celltypes "Endothelial cells,endothelial cell" \
  --tumor_celltypes "Pulmonary vascular smooth muscle cells,smooth muscle cell"

'

####

set -euo pipefail

# this script lives in the root of the module directory
# use this to define other default paths
cd $(dirname "$0")
module_directory=$(pwd)

# define library and project specific arguments
sample_id="SCPCS000490"
library_id="SCPCL000822"
normal_celltypes="Endothelial cells,endothelial cell"
tumor_celltypes="Pulmonary vascular smooth muscle cells,smooth muscle cell"

# build paths to directories for input/output
ref_dir="$module_directory/references"
data_dir="$module_directory/../../data/current/SCPCP000015"
results_dir="$module_directory/results/annotate_tumor_cells_output/$sample_id"

# grab variables from command line
while [ $# -gt 0 ]; do
    if [[ $1 == *'--'* ]]; then
        v="${1/--/}"
        declare $v="$2"
    fi
    shift
done

# define paths to notebooks and scripts run in the workflow
notebook_dir="$module_directory/template_notebooks"
scripts_dir="$module_directory/scripts"
cell_lists_dir="$ref_dir/cell_lists/$sample_id"

# Create table with reference cell types
if [[ -n $normal_celltypes || -n $tumor_celltypes ]]; then

    echo "Saving cell references"

    #Rscript $scripts_dir/select-cell-types.R \
    #  --sce_file "$data_dir/$sample_id/${library_id}_processed.rds" \
    #  --normal_cells "${normal_celltypes}" \
    #  --tumor_cells "${tumor_celltypes}" \
    #  --output_dir "$cell_lists_dir" \
    #  --output_filename "${library_id}_reference-cells.tsv"
#
    reference_cell_file="$cell_lists_dir/${library_id}_reference-cells.tsv"
else
    reference_cell_file=""
fi

echo $module_directory
echo $reference_cell_file

# Step 1: Obtain manual annotations
echo "Getting manual annotations..."
Rscript -e "rmarkdown::render('$notebook_dir/01-marker-genes.Rmd', \
      clean = TRUE, \
      output_dir = '$results_dir', \
      output_file = '${library_id}_marker-gene-report.html', \
      params = list(sample_id = '$sample_id', \
                    library_id = '$library_id', \
                    results_dir = '$results_dir', \
                    reference_cell_file = '$reference_cell_file'), \
      envir = new.env()) \
"
echo "Done"
