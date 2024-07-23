#!/bin/bash

# This script is used to generate the SingleR results for SCPCS000492/SCPCL000824
# The SingleR annotations and tumor cell annotations from SCPCS000490/SCPCL000822 are used as a reference

set -euo pipefail

# Ensure script is being run from its directory
cd $(dirname "$0")
MODULE_DIR="../../.." # navigate up three directories from script to module home

# define input files
DATA_DIR="${MODULE_DIR}/../../data/current/SCPCP000015"
SCE_FILE="${DATA_DIR}/SCPCS000492/SCPCL000824_processed.rds"
REF_SCE_FILE="${DATA_DIR}/SCPCS000490/SCPCL000822_processed.rds"
ANNOTATION_FILE="${MODULE_DIR}/results/annotation_tables/SCPCS000490/SCPCL000822_tumor-classifications.tsv"
OUTPUT_FILE="${MODULE_DIR}/results/annotation_tables/SCPCS000492/SCPCL000824_singler-classifications.tsv"

# run SingleR script
Rscript $MODULE_DIR/scripts/run-singler.R \
  --input_sce_file $SCE_FILE \
  --ref_sce_file $REF_SCE_FILE \
  --ref_annotations_file $ANNOTATION_FILE \
  --output_file $OUTPUT_FILE \
  --threads 4

