#!/bin/bash

# This script is used to set up files needed to run SCimilarity 

# Step 1: The SCimilarity model is downloaded from zenodo (https://zenodo.org/records/10685499) and saved to models/model_v1.1
# Step 2: A reference TSV file is created that includes all possible SCimilarity annotations and associated ontology IDs

# Usage:
#
# ./setup-analysis.sh
#

# By default the full model will be downloaded from zenodo
# To save space, `annotation_only=1` can be used
# This will download only the files needed for CellAnnotation from s3://scpca-references/celltype/scimilarity_references/model_v1.1

# annotation_only=1 ./setup-analysis.sh

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# directory paths
scripts_dir="scripts"
ref_dir="references"

# Define argument defaults
# default is to grab the full model from zenodo, otherwise only the annotation folder is grabbed from S3
annotation_only=${annotation_only:-0}

###################################################################
################## Download SCimilarity model #####################
###################################################################

# output model directory
model_dir="models"
mkdir -p $model_dir

# compressed output file and final model directory
scimilarity_model_dir="${model_dir}/model_v1.1"

if [[ ! -d $scimilarity_model_dir && $annotation_only -eq 0 ]]; then

  # define the url to the model
  scimilarity_model_url="https://zenodo.org/records/10685499/files/model_v1.1.tar.gz?download=1"
  
  # where to store compressed model from zenodo
  compressed_scimilarity_model="${model_dir}/model_v1.1.tar.gz"
  
  echo "Downloading the SCimilarity model..."
  curl -L -o $compressed_scimilarity_model $scimilarity_model_url
  
  echo "Unzipping the model..."
  tar -xzvf $compressed_scimilarity_model -C $scimilarity_model_dir
  rm $compressed_scimilarity_model
  
  echo "Model downloaded and unzipped"
  
elif [[ ! -d $scimilarity_model_dir && $annotation_only -eq 1 ]]; then

  echo "Downloading only the annotation SCimilarity model from s3"
  s3_model='s3://scpca-references/celltype/scimilarity_references/model_v1.1'
  aws s3 cp $s3_model $scimilarity_model_dir --exclude "cellsearch/*" --recursive
  
else
  echo "SCimilarity model already exists"
fi

###################################################################
################## Create cell type reference TSV #################
###################################################################

model_annotations_file="${scimilarity_model_dir}/annotation/reference_labels.tsv"
missing_ontology_file="${ref_dir}/scimilarity-missing-ontology-assignments.tsv"
scimlarity_ontology_file="${ref_dir}/scimilarity-mapped-ontologies.tsv"

Rscript ${scripts_dir}/01-assign-ontology-ids.R \
  --model_annotations_file $model_annotations_file \
  --missing_ontology_tsv $missing_ontology_file \
  --output_ontology_tsv $scimlarity_ontology_file
