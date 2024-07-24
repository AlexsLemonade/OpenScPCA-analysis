#!/usr/bin/env Rscript

# This script generates test data files for use when testing the benchmarking script `../run_doublet-detection-benchmark.sh`.
# Saved in `./data.zip`, the test data files are subsetted version of the datasets used for benchmarking, as described in `../README.md`.
# This script assumes that ../scratch/benchmark-datasets/raw contains the .rds files from <https://doi.org/10.5281/zenodo.4562782>'s `real_datasets.zip` file,
#  and that there is an available `zip` executable on the system.

library(Matrix)
set.seed(2024)

module_base <- rprojroot::find_root(rprojroot::is_renv_project)

input_dir <- file.path(module_base, "scratch", "benchmark-datasets", "raw")
output_dir <- file.path(module_base, "benchmark-test-data", "data")
zip_dir <- glue::glue("{output_dir}.zip")

fs::dir_create(output_dir)

datasets <- c("hm-6k", "pbmc-1B-dm", "pdx-MULTI", "HMEC-orig-MULTI")

# Subset dataset to the following amounts:
n_singlets <- 200
n_doublets <- 25
n_genes <- 500

# Helper function to determine which cells (highest colSum values) or genes to keep (highest rowMean values)
find_names <- function(counts, n, type = "rowmean") {
  if (type == "rowmean") {
    counts <- rowMeans(counts)
  } else if (type == "colsum") {
    counts <- colSums(counts)
  }

  names_to_keep <- counts |>
    sort() |>
    tail(n) |>
    names()

  return(names_to_keep)
}

# Subset and export datasets
datasets |>
  purrr::walk(
    \(dataset) {
      output_file <- file.path(output_dir, glue::glue("{dataset}.rds"))
      raw_data <- readRDS(file.path(input_dir, glue::glue("{dataset}.rds")))

      counts <- raw_data[[1]]
      barcodes <- colnames(counts)
      cell_labels <- raw_data[[2]]
      stopifnot("Different number of labels than barcodes." = length(barcodes) == length(cell_labels))

      # Determine which barcodes to keep - those with highest counts
      keep_barcodes <- c(
        find_names(counts, n_singlets, type = "colsum"),
        find_names(counts, n_doublets, type = "colsum")
      )
      keep_barcode_indices <- which(barcodes %in% keep_barcodes)

      # Determine which genes to keep - those with the highest means
      keep_genes <- find_names(counts, n_genes, type = "rowmean")
      keep_gene_indices <- which(rownames(counts) %in% keep_genes)

      subsetted_data <- list(
        counts[1:n_genes, keep_barcode_indices],
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
