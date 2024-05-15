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
  --sample_metadata sample_metadata.tsv

'
####

set -euo pipefail

# this script lives in the root of the module directory
# use this to define other default paths
cd $(dirname "$0")
module_directory=$(pwd)

# build paths to directories for input/output
workflow_results_dir="$module_directory/results/annotate_tumor_cells_output"
sample_metadata="$module_directory/sample_metadata.tsv"

# define paths to store ref and grab data
ref_dir="$module_directory/references"
data_dir="$module_directory/../../data/current/SCPCP000015"

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

# Run the workflow for each sample in the sample metadata file
# be sure to skip the header
sed 1d $sample_metadata | while read -r sample_id library_id normal_celltypes tumor_celltypes
do

    echo $sample_id

    # define output directory for ref file
    cell_lists_dir="$ref_dir/cell_lists/$sample_id"
    # define output reference file
    reference_cell_file="$cell_lists_dir/${library_id}_reference-cells.tsv"

    # directory to save all results for sample
    sample_results_dir="${workflow_results_dir}/${sample_id}"

    # Create table with reference cell types
    echo "Saving cell references..."
    Rscript $scripts_dir/select-cell-types.R \
      --sce_file "$data_dir/$sample_id/${library_id}_processed.rds" \
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
    echo "Done"
done
