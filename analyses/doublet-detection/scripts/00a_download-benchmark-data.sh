#!/bin/bash

# Download and format ground-truth datasets for use in benchmarking
# Files are saved to ../scratch/benchmark_datasets

set -euo pipefail

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
OUTDIR=../scratch/benchmark_datasets

cd ${SCRIPT_DIR}

# Download and unzip `real_datasets.zip` archive from https://doi.org/10.5281/zenodo.4562782
wget https://zenodo.org/records/4562782/files/real_datasets.zip
unzip real_datasets.zip -d $OUTDIR

# Remove zip file
rm real_datasets.zip

# Run script to read files and export SCE, AnnData files
./00b_format-benchmark-data.R --dir $OUTDIR


