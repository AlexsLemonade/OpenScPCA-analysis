#!/usr/bin/env Rscript

# This script generates test data files for use when testing the benchmarking script `../run_doublet-detection-benchmark.sh`.
# Saved in `./data.zip`, the test data files are subsetted version of the datasets used for benchmarking, as described in `../README.md`.
# This script assumes that ../scratch/benchmark-datasets/raw contains the .rds files from <https://doi.org/10.5281/zenodo.4562782>'s `real_datasets.zip` file and that there is an available `zip` executable on the system.

library(Matrix)
set.seed(2024)

module_base <- rprojroot::find_root(rprojroot::is_renv_project)

input_dir <- file.path(module_base, "scratch", "benchmark-datasets", "raw")
output_dir <- file.path(module_base, "benchmark-test-data", "data")
zip_dir <- glue::glue("{output_dir}.zip")

fs::dir_create(output_dir)

datasets <- c("hm-6k", "pbmc-1B-dm", "pdx-MULTI", "HMEC-orig-MULTI")

# For each dataset, subset to 50 droplets, including 45 singlets and 5 doublets, and the first 100 genes
# In all original datasets, doublets have been pre-sorted to the end of the list
datasets |>
  purrr::walk(
    \(dataset) {
      output_file <- file.path(output_dir, glue::glue("{dataset}_subset.rds"))

      raw_data <- readRDS(file.path(input_dir, glue::glue("{dataset}.rds")))

      # keep 45 random singlets and 5 random doublets
      counts <- raw_data[[1]]
      barcodes <- colnames(counts)
      cell_labels <- raw_data[[2]]
      stopifnot("Different number of labels than barcodes." = length(barcodes) == length(cell_labels))

      keep_barcodes <- c(
        sample(barcodes[cell_labels == "singlet"], 45),
        sample(barcodes[cell_labels == "doublet"], 5)
      )
      keep_barcode_indices <- which(barcodes %in% keep_barcodes)

      subsetted_data <- list(
        counts[1:100, keep_barcode_indices],
        cell_labels[keep_barcode_indices]
      )

      saveRDS(subsetted_data, output_file)
    }
  )

# Compress the output directory into a zip file and remove the output directory
zip(
  zipfile = zip_dir,
  files = dir(output_dir, full.names = TRUE),
  flags = "-rj"
)
fs::dir_delete(output_dir)
