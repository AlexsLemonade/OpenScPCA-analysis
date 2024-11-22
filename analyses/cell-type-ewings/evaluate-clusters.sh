#!/bin/bash

# This script is used to evaluate clustering with varying parameters for all samples in SCPCP000015
# Clusters are calculated using the following parameters: 

# Louvain: 
#  - nearest neighbors: 5, 10, 15, 20, 25, 30, 35, 40
#  - resolution: 0.5, 1, 1.5

# Leiden CPM: 
#  - nearest neighbors: 5, 10, 15, 20, 25, 30, 35, 40
#  - resolution: 0.001, 0.005, 0.01

# Leiden modularity:
#  - nearest neighbors: 5, 10, 15, 20, 25, 30, 35, 40
#  - resolution: 0.5, 1, 1.5

# For each sample, all clustering results are saved as a TSV file and metrics are summarized in a summary HTML report 

# Usage: ./evaluate-clusters.sh

set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")
module_dir=$(pwd)

# use 4 CPUs and define seed
threads=4
seed=2024

# set up input and output file paths
# make sure full paths are provided for notebook rendering
data_dir="${module_dir}/../../data/current/SCPCP000015"
results_dir="${module_dir}/results/clustering"
mkdir -p ${results_dir}

# define scripts and notebook directories
scripts_dir="scripts/clustering-workflow"
notebook_dir="template_notebooks/clustering-workflow"

# define all sample IDs 
sample_ids=$(basename -a ${data_dir}/SCPCS*)
# perform clustering on each library 
for sample_id in $sample_ids; do

    # make sure sample results directory exists
    sample_results_dir="${results_dir}/${sample_id}"
    mkdir -p $sample_results_dir

    for sce_file in $data_dir/$sample_id/*_processed.rds; do

    # define library ID
        library_id=$(basename $sce_file | sed 's/_processed.rds$//')

        # define output files
        cluster_results="${sample_results_dir}/${library_id}_cluster-results.tsv"

        # generate clusters 
        # use all three clustering methods
        echo "Generating clusters for sample: $sample_id and library: $library_id"
        Rscript $scripts_dir/01-clustering.R \
            --sce_file $sce_file \
            --output_file $cluster_results \
            --louvain --leiden_mod --leiden_cpm \
            --threads $threads \
            --seed $seed

        # render notebook
        Rscript -e "rmarkdown::render('$notebook_dir/01-clustering-metrics.Rmd', \
          clean = TRUE, \
          output_dir = '$sample_results_dir', \
          output_file = '${library_id}_cluster-summary-report.html', \
          params = list(library_id = '$library_id', \
                        sce_file = '$sce_file', \
                        cluster_results_file = '$cluster_results'), \
          envir = new.env()) \
        "
    done
done
