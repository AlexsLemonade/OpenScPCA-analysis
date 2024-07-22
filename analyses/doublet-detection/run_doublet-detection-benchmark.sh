#!/bin/bash

# This script runs the benchmarking portion of the `doublet-detection` module
# Usage: ./run_doublet-detection-benchmark.sh
# Additional arguments:
# The number of cores (default is 4) to use with scDblFinder can be set with the `--cores` option:
#    cores=2 ./run_doublet-detection-benchmark.sh
# If this script is being run in test mode, use `--test=1` option:
#    test=1 ./run_doublet-detection-benchmark.sh

set -euo pipefail

# input variables
CORES=${cores:-4}
TEST=${test:-0}


# Ensure script is being run from its directory
MODULE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd ${MODULE_DIR}

TEMPLATE_NB_DIR="template-notebooks" # directory with template notebooks
EXPLORE_NB_DIR="exploratory-notebooks"  # directory with exploratory notebooks

# Define and create benchmark directories
DATA_DIR="scratch/benchmark-datasets"
RESULTS_DIR="results/benchmark-results"
RESULTS_NB_DIR="${RESULTS_DIR}/rendered-notebooks"
mkdir -p ${DATA_DIR}
mkdir -p ${RESULTS_DIR}
mkdir -p ${RESULTS_NB_DIR}

# Remove any existing data and results to ensure there are no conflicts
rm -rf ${DATA_DIR}/*
rm -rf ${RESULTS_DIR}/*
rm -rf ${RESULTS_NB_DIR}/*

# define benchmarking dataset names
bench_datasets=("hm-6k" "pbmc-1B-dm" "pdx-MULTI" "HMEC-orig-MULTI")

# download benchmarking files to $DATA_DIR/raw
if [[ $TEST -eq 1 ]]; then
    # Download the subsetted test datasets
    curl -O https://github.com/AlexsLemonade/OpenScPCA-analysis/raw/c9aad22478c3888089bcff9d6a04f011844c5256/analyses/doublet-detection/benchmark-test-data/data.zip
    unzip data.zip -d ${DATA_DIR}/raw
    rm data.zip
else
    # Download and unzip `real_datasets.zip` archive from https://doi.org/10.5281/zenodo.4562782
    curl -O https://zenodo.org/records/4562782/files/real_datasets.zip
    unzip real_datasets.zip -d ${DATA_DIR}/raw
    rm real_datasets.zip
fi

for dataset in "${bench_datasets[@]}"; do

    echo "============================${dataset}=========================="

    # formatted SCE and AnnData files will be saved here
    DATASET_DIR=${DATA_DIR}/$dataset
    mkdir -p $DATASET_DIR

    # Read raw downloaded data and export SCE, AnnData files
    ./scripts/00_format-benchmark-data.R --dataset ${dataset} --input_dir ${DATA_DIR}/raw --output_dir ${DATASET_DIR}

    # Infer doublets with scDblFinder
    SCDBLFINDER_TSV=${RESULTS_DIR}/${dataset}_scdblfinder.tsv
    ./scripts/01a_run-scdblfinder.R --input_sce_file ${DATASET_DIR}/${dataset}.rds --output_file ${SCDBLFINDER_TSV} --cores $CORES --benchmark

    # Infer doublets with scrublet
    SCRUBLET_TSV=${RESULTS_DIR}/${dataset}_scrublet.tsv
    ./scripts/01b_run-scrublet.py --input_anndata_file ${DATASET_DIR}/${dataset}.h5ad --output_file ${SCRUBLET_TSV}

    # Explore each individual set of doublet results
    Rscript -e "rmarkdown::render('${TEMPLATE_NB_DIR}/02_explore-benchmark-results.Rmd',
            output_dir = '${RESULTS_NB_DIR}',
            output_file = '${dataset}_doublet-results.html',
            params = list(dataset = '${dataset}'),
            clean = TRUE)"
done

# Compare doublet inferences across methods, on all datasets processed
Rscript -e "rmarkdown::render('${EXPLORE_NB_DIR}/03_compare-benchmark-results.Rmd',
            output_dir = '${RESULTS_NB_DIR}',
            output_file = 'compare-benchmark-results.nb.html',
            clean = TRUE)"
