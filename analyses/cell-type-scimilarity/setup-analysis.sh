#!/bin/bash

# This script is used to set up files needed to run SCimilarity 

# Step 1: The SCimilarity model is downloaded from zenodo (https://zenodo.org/records/10685499) and saved to models/model_v1.1
# Step 2: A reference TSV file is created that includes all possible SCimilarity annotations and associated ontology IDs

# Usage:
#
# ./run-analysis.sh
#

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# directory paths
scripts_dir="scripts"
ref_dir="references"

###################################################################
################## Download SCimilarity model #####################
###################################################################

# define the url to the model
scimilarity_model_url="https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1"

# output model directory
model_dir="models"
mkdir -p $model_dir

# compressed output file and final model directory
compressed_scimilarity_model="${model_dir}/model_v1.1.tar.gz"
scimilarity_model_dir="${model_dir}/model_v1.1"

if [[ ! -d $scimilarity_model_dir ]]; then
  echo "Downloading the SCimilarity model..."
  curl -L -o $compressed_scimilarity_model $scimilarity_model_url
  echo "Unzipping the model..."
  tar -xzvf $compressed_scimilarity_model -C $scimilarity_model_dir
  rm $compressed_scimilarity_model
  echo "Model downloaded and unzipped"
else
  echo "SCimilarity model already exists"
fi

###################################################################
################## Create cell type reference TSV #################
###################################################################

missing_ontology_file="${ref_dir}/scimilarity-missing-ontology-assignments.tsv"
scimlarity_ontology_file="${ref_dir}/scimilarity-mapped-ontologies.tsv"

Rscript ${scripts_dir}/01-assign-ontology-ids.R \
  --model_dir $scimilarity_model_dir \
  --missing_ontology_tsv $missing_ontology_file \
  --output_ontology_tsv $scimlarity_ontology_file
