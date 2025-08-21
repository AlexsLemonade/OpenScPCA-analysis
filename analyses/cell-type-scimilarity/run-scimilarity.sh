#!/bin/bash

# This script is used to annotate cell types with SCimilarity for all libraries in a single ScPCA project 
# A single TSV file will be saved for every library with the annotations and min_dist reported by SCimilarity 

# Prior to running this scrip, the model must be downloaded directly the setup-analysis.sh script

# Usage: ./run-scimilarity.sh SCPCP000001

set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")

project_id=$1

data_dir="../../data/current"

# directory to store scimilarity results 
scimilarity_results_dir="results/scimilarity/${project_id}" 
mkdir -p ${scimilarity_results_dir}

# define path to input model 
model_dir="models/model_v1.1"

for sample_dir in ${data_dir}/${project_id}/SCPCS*; do

  # grab sample id
  sample_id=$(basename $sample_dir)

  # define output folder and make sure it exists
  sample_results_dir="${scimilarity_results_dir}/${sample_id}"
  mkdir -p ${sample_results_dir}

  # run script to export tsv file on all processed objects
  for anndata_file in $sample_dir/*_processed_rna.h5ad; do

    # define library ID
    library_id=$(basename ${anndata_file%_processed_rna.h5ad})

    # define output file
    scimilarity_predictions_file="${sample_results_dir}/${library_id}_scimilarity-predictions.tsv.gz"
    
    echo "Assigning cell types for ${library_id}"
    python3 scripts/02-run-scimilarity.py \
      --model_dir $model_dir \
      --processed_h5ad_file $anndata_file \
      --predictions_tsv $scimilarity_predictions_file 

  done 

done
