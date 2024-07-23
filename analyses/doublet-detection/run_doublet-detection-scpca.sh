#!/bin/bash

# This script runs doublet detection over ScPCA data for a given project
# Usage: ./run_doublet-detection-scpca.sh {scpca project id}
# Additional arguments:
# The number of cores (default is 4) to use with scDblFinder can be set with the `--cores` option:
#    cores=2 ./run_doublet-detection-benchmark.sh {scpca project id}

set -euo pipefail

# input variables
CORES=${cores:-4}

# Ensure script is being run from its directory
MODULE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd ${MODULE_DIR}

PROJECT_ID=$1

# Define directories
DATA_DIR="../../data/current"
RESULTS_DIR="results/scpca-results/${PROJECT_ID}"
mkdir -p ${RESULTS_DIR}

# Detect doublets on each processed SCE file in each sample directory
for SAMPLE_DIR in ${DATA_DIR}/${PROJECT_ID}/SCPCS*; do
    SAMPLE_ID=$(basename $SAMPLE_DIR)
    echo "Processing ${SAMPLE_ID}..."

    SAMPLE_RESULTS_DIR=${RESULTS_DIR}/${SAMPLE_ID}
    mkdir -p ${SAMPLE_RESULTS_DIR}

    for SCE_FILE in ${SAMPLE_DIR}/*_processed.rds; do
        TSV_FILE=$(basename "${SCE_FILE%.rds}_scdblfinder.tsv")
        Rscript scripts/01a_run-scdblfinder.R \
            --input_sce_file ${SCE_FILE} \
            --output_file ${SAMPLE_RESULTS_DIR}/${TSV_FILE} \
            --cores $CORES
    done
done
