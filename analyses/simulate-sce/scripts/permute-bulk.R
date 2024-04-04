#!/usr/bin/env Rscript

# Script to generate permuted bulk RNAseq data from a real dataset

# Works on a single bulk file, and permutes the values within each gene


# Parse arguments --------------------------------------------------------------

library(optparse)

# Define and parse the command line arguments
option_list <- list(
  make_option(
    c("-f", "--bulk_file"),
    type = "character",
    default = NULL,
    help = "The sample directory for the real data to use for simulation"
  ),
  make_option(
    c("-o", "--output_dir"),
    type = "character",
    default = NULL,
    help = "The output directory. Output files will be given the same names as input"
  ),
  make_option(
    c("--seed"),
    type = "integer",
    default = 2024,
    help = "Random number seed for simulation."
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# set seed
set.seed(opts$seed)

# check outputs
stopifnot(
  "Could not find `bulk_file`" = file.exists(opts$bulk_file),
  "Could not find `output_dir`" = dir.exists(opts$output_dir)
)

# read in data table & convert to matrix
bulk_data <- readr::read_tsv(
  opts$bulk_file,
  col_types = readr::cols(gene_id = "c", .default = "d")
) |>
  tibble::column_to_rownames(var = "gene_id") |>
  as.matrix()

stopifnot("Data table must have at least 2 samples" = ncol(bulk_data) >= 2)

# permute by row (comes out transposed, so transpose back and add sample names)
permuted_data <- apply(bulk_data, 1, sample) |>
  t()
colnames(permuted_data) <- colnames(bulk_data)

# convert back to data frame
permuted_data <- tibble::as_tibble(permuted_data, rownames = "gene_id")

# write out permuted data

output_file <- file.path(opts$output_dir, basename(opts$bulk_file))
readr::write_tsv(permuted_data, file = output_file)
