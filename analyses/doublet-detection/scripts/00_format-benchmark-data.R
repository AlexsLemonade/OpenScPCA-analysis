#!/usr/bin/env Rscript

# This script converts the datasets used for benchmarking into SCE and AnnData objects,
# mirroring ScPCA dataset format.
# The original data files were obtained from Zenodo and are named <dataset name>.rds.
# Each is an RDS object containing a list of two items:
#  [[1]] A raw counts matrix
#  [[2]] A vector of doublet/singlet calls for each barcode
# The exported files are named `<dataset name>-<sce/anndata>.<rds/h5ad>.
# Each has a variable `ground_truth_doublets` representing the singlet/doublet calls:
#   - In the SCE file, this is in the colData slot
#   - In the AnnData file, this is in the obs slot

# Load renv
project_root <- here::here()
renv::load(project_root)


library(optparse)
library(SingleCellExperiment)
library(zellkonverter)


option_list <- list(
  make_option(
    c("-d", "--dir"),
    type = "character",
    help = "Directory containing data files to format."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))



# We're only focusing on these datasets for benchmarking
datanames <- c("hm-6k", "pbmc-1B-dm", "pdx-MULTI", "HMEC-orig-MULTI")

# Perform conversion
datanames |>
  purrr::walk(
    \(dataname) {
      input_file <- file.path(opts$dir, glue::glue("{dataname}.rds"))
      output_sce_file <- file.path(opts$dir, glue::glue("{dataname}_sce.rds"))
      output_anndata_file <- file.path(opts$dir, glue::glue("{dataname}_anndata.h5ad"))

      dat <- readRDS(input_file)
      mat <- dat[[1]]
      calls <- dat[[2]]

      # Create and export SCE
      sce <- SingleCellExperiment(assays = list(counts = mat))
      sce$ground_truth_doublets <- calls
      readr::write_rds(sce, output_sce_file)

      # Export AnnData version
      writeH5AD(sce, output_anndata_file)
    }
  )
