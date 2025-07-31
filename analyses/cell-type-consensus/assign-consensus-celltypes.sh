#!/bin/bash

# This script is used to assign consensus cell types to all cells in a given ScPCA project  
# The existing cell type annotations from SingleR and CellAssign are saved alongside the consensus cell types
# A single TSV file will be saved for every library 

# Usage: ./assign-consensus-celltypes.sh SCPCP000001


set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")

project_id=$1

data_dir="../../data/current"

# directory to store consensus results 
consensus_results_dir="results/cell-type-consensus/${project_id}" 
mkdir -p ${consensus_results_dir}

# define reference input files 
blueprint_ref_file="references/blueprint-mapped-ontologies.tsv"
panglao_ref_file="references/panglao-cell-type-ontologies.tsv"
consensus_ref_file="references/consensus-cell-type-reference.tsv"

# marker gene file to use for validation 
validation_markers_file="references/validation-markers.tsv"

for sample_dir in ${data_dir}/${project_id}/SCPCS*; do

  # grab sample id
  sample_id=$(basename $sample_dir)

  # define output folder and make sure it exists
  sample_results_dir="${consensus_results_dir}/${sample_id}"
  mkdir -p ${sample_results_dir}

  # run script to export tsv file on all processed objects
  for sce_file in $sample_dir/*_processed.rds; do

    # define library ID
    library_id=$(basename $sce_file | sed 's/_processed.rds$//')

    # define output files
    consensus_output_file="${sample_results_dir}/${library_id}_consensus-cell-type-assignments.tsv.gz"
    gene_exp_output_file="${sample_results_dir}/${library_id}_marker-gene-expression.tsv.gz"
    
    echo "Assigning cell types for ${library_id}"
    Rscript scripts/04-assign-consensus-celltypes.R \
      --sce_file $sce_file \
      --blueprint_ref_file $blueprint_ref_file\
      --panglao_ref_file $panglao_ref_file \
      --consensus_ref_file $consensus_ref_file \
      --marker_gene_file $validation_markers_file \
      --consensus_output_file $consensus_output_file \
      --gene_exp_output_file $gene_exp_output_file

  done 

done
