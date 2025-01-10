#!/bin/bash

# This script is used to create a single table with cell type assignments all cells from all ScPCA samples  
# The existing cell type annotations from SingleR and CellAssign are saved to a TSV file for each sample 
# Then all TSV files are combined into a single file and consensus cell types are assigned 

# Usage: ./assign-consensus-celltypes.sh 


set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")
#module_dir=$(pwd)

data_dir="../../data/current"
# path to save consensus results 
scpca_consensus_assignments_file="results/scpca-consensus-celltype-assignments.tsv.gz"
# directory to store all individual tsv files 
celltype_tsv_dir="results/original-celltype-assignments" 
mkdir -p ${celltype_tsv_dir}

# define reference input files 
panglao_ref_file="references/panglao-cell-type-ontologies.tsv"
consensus_ref_file="references/consensus-cell-type-reference.tsv"

# run script to export tsv file on all processed objects
for sce_file in $data_dir/SCPCP*/SCPCS*/*_processed.rds; do

    # define library ID
    library_id=$(basename $sce_file | sed 's/_processed.rds$//')
    
    echo "Grabbing cell types for ${library_id}"
    # get celltypes as tsv file 
    Rscript scripts/03-save-coldata.R \
      --sce_file $sce_file \
      --output_file ${celltype_tsv_dir}/${library_id}_celltype-assignments.tsv

done 

echo "Combining TSVs and adding consensus labels"
# run script to combine all tsv files and assign consensus cell types
Rscript scripts/04-combine-celltype-tables.R \
  --celltype_tsv_dir $celltype_tsv_dir \
  --panglao_ref_file $panglao_ref_file \
  --consensus_ref_file $consensus_ref_file \
  --output_file $scpca_consensus_assignments_file
