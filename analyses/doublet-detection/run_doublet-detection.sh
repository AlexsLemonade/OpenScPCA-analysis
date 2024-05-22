#!/bin/bash

set -euo pipefail

# Set up --------------

# Ensure script is being run from its directory
MODULE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd ${MODULE_DIR}

# Create directories
mkdir -p scratch/benchmark_datasets
mkdir -p results/benchmark_results


# Step 1: Benchmark doublet detection methods on ground truth datasets -----------

# Download and format data
./scripts/00a_download-benchmark-data.sh

# Infer doublets with scDblFinder
./scripts/01a_detect-doublets.R

# Infer doublets with scrublet 01b_detect-doublets.py
./scripts/01b_detect-doublets.py
