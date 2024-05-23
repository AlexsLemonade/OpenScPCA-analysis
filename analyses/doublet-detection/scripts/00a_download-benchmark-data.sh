#!/bin/bash

# Download and format ground-truth datasets for use in benchmarking from https://doi.org/10.5281/zenodo.4562782 (`real_datasets`)
# This script takes one argument, the output directory where data should be saved

set -euo pipefail

SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
OUTDIR=$1

cd ${SCRIPT_DIR}

# Download and unzip `real_datasets.zip` archive from https://doi.org/10.5281/zenodo.4562782
wget https://zenodo.org/records/4562782/files/real_datasets.zip
unzip real_datasets.zip -d $OUTDIR

# Remove zip file
rm real_datasets.zip

# Run script to read files and export SCE, AnnData files
./00b_format-benchmark-data.R --dir $OUTDIR


