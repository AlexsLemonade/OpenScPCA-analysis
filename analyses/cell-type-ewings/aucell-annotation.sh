#!/bin/bash

# This script identfies tumor cells using AUCell for all samples in SCPCP000015
# Usage: ./aucell-annotation.sh

set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")
module_dir=$(pwd)

# use 4 CPUs
threads=4

# set up input and output file paths
# make sure full paths are provided for notebook rendering
data_dir="${module_dir}/../../data/current/SCPCP000015"
results_dir="${module_dir}/results/aucell_annotation"
mkdir -p ${results_dir}

# define path to marker genes file
marker_genes_file="${module_dir}/references/tumor-marker-genes.tsv"

# define scripts and notebook directories
scripts_dir="scripts/auc-workflow"
notebook_dir="template_notebooks/auc-workflow"

# define reference sample inputs and outputs
ref_sample="SCPCS000490"
ref_library="SCPCL000822"
ref_sce="${data_dir}/${ref_sample}/${ref_library}_processed.rds"

ref_results_dir="${results_dir}/${ref_sample}"
mkdir -p $ref_results_dir

ref_auc_results="${results_dir}/${ref_sample}/${ref_library}_auc-classifications.tsv"

echo "Running AUCell for reference sample: $ref_sample and library: $ref_library"
# first run AUCell on reference sample
auc_threshold=$(Rscript $scripts_dir/01-run-aucell.R \
    --sce_file $ref_sce \
    --output_file $ref_auc_results \
    --return_auc)

# now run AUCell, gene set scores, and generate AUCell report for each sample that is not the ref sample
sample_ids=$(basename -a ${data_dir}/SCPCS*)
for sample_id in $sample_ids; do

    # make sure sample results directory exists
    sample_results_dir="${results_dir}/${sample_id}"
    mkdir -p $sample_results_dir

    for sce_file in $data_dir/$sample_id/*_processed.rds; do

        # define library ID
        library_id=$(basename $sce_file | sed 's/_processed.rds$//')

        # define output files
        auc_results="${sample_results_dir}/${library_id}_auc-classifications.tsv"
        geneset_results="${sample_results_dir}/${library_id}_gene-set-scores.tsv"
        marker_gene_results=${sample_results_dir}/${library_id}_marker-gene-classifications.tsv"

        # only run AUCell if sample is NOT the ref sample
        if [[ ($sample_id != $ref_sample) && ($library_id != $ref_library) && (! -f $auc_results) ]]; then

            echo "Running AUCell for sample: $sample_id and library: $library_id"

            # run AUCell
            Rscript $scripts_dir/01-run-aucell.R \
                --sce_file $sce_file \
                --auc_threshold $auc_threshold \
                --output_file $auc_results \
                --threads $threads

        fi

        # generate gene set scores
        Rscript $scripts_dir/02-calculate-gene-set-scores.R \
          --sce_file $sce_file \
          --output_file $geneset_results

        # render notebook
        Rscript -e "rmarkdown::render('$notebook_dir/01-auc-results.Rmd', \
          clean = TRUE, \
          output_dir = '$sample_results_dir', \
          output_file = '${library_id}_aucell-report.html', \
          params = list(library_id = '$library_id', \
                        sce_file = '$sce_file', \
                        aucell_results_file = '$auc_results', \
                        ref_auc_results = '$ref_auc_results', \
                        geneset_scores_file = '$geneset_results', \
                        marker_genes_file = '$marker_genes_file', \
                        output_file = '${marker_gene_results}', \
                        auc_threshold = $auc_threshold), \
          envir = new.env()) \
        "
    done
done
