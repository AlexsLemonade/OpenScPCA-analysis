#!/usr/bin/env Rscript

# This script generates a test data file for use when testing the benchmarking script `../run_doublet-detection-benchmark.sh`.
# The test data file is a subsetted version of the `hm-6k` dataset used for benchmarking, as described in `../README.md`.
# This script assumes that ../scratch/benchmark-datasets/raw contains benchmarking datasets <https://doi.org/10.5281/zenodo.4562782>


library(Matrix)

input_file <- file.path("..", "scratch", "benchmark-datasets", "raw", "hm-6k.rds")
output_file <- file.path("hm-6k_subset.rds")


raw_data <- readRDS(input_file)

# Keep 50 droplets, including both singlets and doublets, and 100 genes
subsetted_data <- list(
    raw_data[[1]][c(1:45, 6802:6806), 1:100],
    raw_data[[2]][c(1:45, 6802:6806)]
)

saveRDS(subsetted_data, output_file)
