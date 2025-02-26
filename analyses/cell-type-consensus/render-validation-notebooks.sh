#!/bin/bash

# This script is used to create an html notebook summarizing the expression of cell type specific marker genes in the assigned consensus cell types 
# One notebook will pe produced for each project
# The input consensus cell types are obtained from running the cell-type-consensus module in OpenScPCA-nf 
# and can be downloaded using the download script in the root directory of the repo using the following command: 

# ./download-results.py --module cell-type-consensus 

# To use this script run:

#  ./render-validation-notebooks.sh


set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")

# directory with results from cell-type-consensus
consensus_results_dir="../../data/current/results/cell-type-consensus"

# directory to store rendered notebooks 
output_notebook_dir="exploratory-notebooks/cell-type-validation-notebooks" 
mkdir -p ${output_notebook_dir}

# render one notebook for each project 
for project_id in $(basename ${consensus_results_dir}/SCPCP*); do

    output_file="${output_notebook_dir}/${project_id}_cell-type-validation-summary.html"

    # skip the large project for now 
    if [[ (! -f $output_file) && ($project_id != "SCPCP000003") && ($project_id != "SCPCP000008")]]; then

        echo "Rendering notebook for $project_id"
        Rscript -e "rmarkdown::render('template-notebooks/marker-gene-validation.Rmd', \
          clean = TRUE, \
          output_dir = '$output_notebook_dir', \
          output_file = '${project_id}_cell-type-validation-summary.html', \
          params = list(project_id = '$project_id'), \
          envir = new.env()) \
        "
    fi

done
