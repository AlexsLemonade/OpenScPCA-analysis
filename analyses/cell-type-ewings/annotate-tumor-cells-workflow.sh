#!/usr/bin/env bash

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
