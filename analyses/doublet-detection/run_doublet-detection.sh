#!/bin/bash

set -euo pipefail

# Set up --------------

# Ensure script is being run from its directory
MODULE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd ${MODULE_DIR}

TEMPLATE_NB_DIR=template-notebooks # directory with template notebooks
EXPLORE_NB_DIR=exploratory-notebooks  # directory with exploratory notebooks
###########################################################################################
########## Step 1: Benchmark doublet detection methods on ground truth datasets ###########
###########################################################################################

# Create benchmark directories
BENCH_DATA_DIR=scratch/benchmark-datasets
BENCH_RESULTS_DIR=results/benchmark-results
BENCH_TEMPLATE_NB_DIR=${BENCH_RESULTS_DIR}/rendered-notebooks
mkdir -p ${BENCH_DATA_DIR}
mkdir -p ${BENCH_RESULTS_DIR}
mkdir -p ${BENCH_TEMPLATE_NB_DIR}


# define benchmarking datasets to use
bench_datasets=("hm-6k" "pbmc-1B-dm" "pdx-MULTI" "HMEC-orig-MULTI")

# Download and unzip `real_datasets.zip` archive from https://doi.org/10.5281/zenodo.4562782
# Files are saved in $BENCH_DATA_DIR/raw
wget https://zenodo.org/records/4562782/files/real_datasets.zip
unzip real_datasets.zip -d ${BENCH_DATA_DIR}/raw
rm real_datasets.zip

for dataset in "${bench_datasets[@]}"; do

    # formatted SCE and AnnData files will be saved here
    DATASET_DIR=${BENCH_DATA_DIR}/$dataset
    mkdir -p $DATASET_DIR

    # Read raw downloaded data and export SCE, AnnData files
    ./scripts/00_format-benchmark-data.R --dataset_name ${dataset} --input_dir ${BENCH_DATA_DIR}/raw --output_dir ${DATASET_DIR}

    # Infer doublets with scDblFinder
    ./scripts/01a_run-scdblfinder.R --dataset_name ${dataset} --data_dir ${DATASET_DIR} --results_dir ${BENCH_RESULTS_DIR}

    # Infer doublets with scrublet
    ./scripts/01b_run-scrublet.py --dataset_name ${dataset} --data_dir ${DATASET_DIR} --results_dir ${BENCH_RESULTS_DIR}

    # Explore each individual set of doublet results
    Rscript -e "rmarkdown::render('${TEMPLATE_NB_DIR}/02_explore-benchmark-results.Rmd',
            output_dir = '${BENCH_TEMPLATE_NB_DIR}',
            output_file = '${dataset}-doublet-results.html',
            params = list(dataset = '${dataset}'),
            clean = TRUE)"
done

# Compare doublet inferences across methods, on all datasets processed
Rscript -e "rmarkdown::render('${EXPLORE_NB_DIR}/03_compare-benchmark-results.Rmd', clean = TRUE)"
