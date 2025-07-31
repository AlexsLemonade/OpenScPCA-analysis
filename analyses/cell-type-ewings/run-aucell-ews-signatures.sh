#!/bin/bash

# This script is used to get AUC values for a set of EWS-FLI1 high/low gene signatures using AUCell
# The input gene signatures are those saved in `references/gene_signatures` and MSigDB gene sets noted in `references/msigdb-gene-sets.tsv`

# AUCell is run on each library and the merged object 
# The output is a TSV for each library/object with the following columns: 
# `gene_set`, `barcodes`, `auc`, and `auc_threshold`

# Usage: ./run-aucell-ews-signatures.sh

# By default AUCell is run with an aucMaxRank of 1% (max_rank_threshold=0.01) of the total genes detected in a sample
# To change this use the following command, specifying a max_rank_threshold between 0-1

# max_rank_threshold=.05 ./run-aucell-ews-signatures.sh

set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")
module_dir=$(pwd)

# use 4 CPUs and define seed
threads=${threads:-4}
seed=${seed:-2025}

# set auc max rank value
max_rank_threshold=${max_rank_threshold:-0.01}

# set up input and output file paths
data_dir="../../data/current/SCPCP000015"
merged_sce_file="../../data/current/results/merge-sce/SCPCP000015/SCPCP000015_merged.rds"

results_dir="results/aucell-ews-signatures"
mkdir -p ${results_dir}

merged_results_file="${results_dir}/SCPCP000015_auc-ews-gene-signatures.tsv"

# gene signature directory 
gene_signatures_dir="${module_dir}/references/gene_signatures"
msigdb_geneset_file="${module_dir}/references/msigdb-gene-sets.tsv"

# define all sample IDs 
sample_ids=$(basename -a ${data_dir}/SCPCS*)
# run AUCell on each library 
for sample_id in $sample_ids; do

    # make sure sample results directory exists
    sample_results_dir="${results_dir}/${sample_id}"
    mkdir -p $sample_results_dir

    for sce_file in $data_dir/$sample_id/*_processed.rds; do

        # define library ID
        library_id=$(basename $sce_file | sed 's/_processed.rds$//')

        # define output file 
        auc_results_file="${sample_results_dir}/${library_id}_auc-ews-gene-signatures.tsv"

        # run AUCell
        echo "Running AUCell for $library_id"
        Rscript scripts/aucell-ews-signatures/01-aucell.R \
          --sce_file $sce_file \
          --custom_geneset_dir $gene_signatures_dir \
          --msigdb_genesets $msigdb_geneset_file \
          --max_rank_threshold $max_rank_threshold \
          --output_file $auc_results_file \
          --threads $threads \
          --seed $seed
      done
done


# run AUCell on merged object 
echo "Running AUCell for merged object" 
Rscript scripts/aucell-ews-signatures/01-aucell.R \
    --sce_file $merged_sce_file \
    --custom_geneset_dir $gene_signatures_dir \
    --msigdb_genesets $msigdb_geneset_file \
    --max_rank_threshold $max_rank_threshold \
    --is_merged \
    --output_file $merged_results_file \
    --threads $threads \
    --seed $seed
