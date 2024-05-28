#!/bin/bash

set -euo pipefail

# Set up --------------

# Ensure script is being run from its directory
MODULE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd ${MODULE_DIR}

###########################################################################################
########## Step 1: Benchmark doublet detection methods on ground truth datasets ###########
###########################################################################################

# Create benchmark directories
BENCH_DATA_DIR=scratch/benchmark_datasets
BENCH_RESULTS_DIR=results/benchmark_results
mkdir -p ${BENCH_DATA_DIR}
mkdir -p ${BENCH_RESULTS_DIR}

# define benchmarking datasets to use
bench_datasets=("hm-6k" "pbmc-1B-dm" "pdx-MULTI" "HMEC-orig-MULTI")

# Download and unzip `real_datasets.zip` archive from https://doi.org/10.5281/zenodo.4562782
# Files are saved in $BENCH_DATA_DIR/raw
#wget https://zenodo.org/records/4562782/files/real_datasets.zip
#unzip real_datasets.zip -d ${BENCH_DATA_DIR}/raw
#rm real_datasets.zip

for dataset in "${bench_datasets[@]}"; do

    # formatted SCE and AnnData files will be saved here
    DATASET_DIR=${BENCH_DATA_DIR}/$dataset
    mkdir -p $DATASET_DIR

    # Read raw downloaded data and export SCE, AnnData files
    ./scripts/00_format-benchmark-data.R --dataset_name ${dataset} --input_dir ${BENCH_DATA_DIR}/raw --output_dir ${DATASET_DIR}

    # Infer doublets with scDblFinder
    #./scripts/01a_run-scdblfinder.R --dataset_name ${dataset} --data_dir ${DATASET_DIR} --results_dir ${BENCH_RESULTS_DIR}

    # Infer doublets with scrublet
   # ./scripts/01b_run-scrublet.py --dataset_name ${dataset} --data_dir ${DATASET_DIR} --results_dir ${BENCH_RESULTS_DIR}

done
