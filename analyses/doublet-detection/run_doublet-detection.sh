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
BENCH_DATA_DIR_RAW=scratch/benchmark_datasets
BENCH_DATA_DIR=${BENCH_DATA_DIR_RAW}/formatted
BENCH_RESULT_DIR=results/benchmark_results
mkdir -p ${BENCH_DATA_DIR} # also creates ${BENCH_DATA_DIR_RAW}
mkdir -p ${BENCH_RESULT_DIR}

# define benchmarking datasets to use
bench_datasets=("hm-6k" "pbmc-1B-dm" "pdx-MULTI" "HMEC-orig-MULTI")

# Download and unzip `real_datasets.zip` archive from https://doi.org/10.5281/zenodo.4562782
wget https://zenodo.org/records/4562782/files/real_datasets.zip
unzip real_datasets.zip -d ${BENCH_DATA_DIR_RAW}
rm real_datasets.zip

for dataset in "${bench_datasets[@]}"; do

    # Read raw downloaded data and export SCE, AnnData files
    ./scripts/00_format-benchmark-data.R --dataset_name ${dataset} --input_dir ${BENCH_DATA_DIR_RAW} --output_dir ${BENCH_DATA_DIR}

    # Infer doublets with scDblFinder
    ./scripts/01a_detect-doublets.R --dataset_name ${dataset} --data_dir ${BENCH_DATA_DIR} --results_dir ${BENCH_RESULT_DIR}

    # Infer doublets with scrublet
    ./scripts/01b_detect-doublets.py --dataset_name ${dataset} --data_dir ${BENCH_DATA_DIR} --results_dir ${BENCH_RESULT_DIR}

done
