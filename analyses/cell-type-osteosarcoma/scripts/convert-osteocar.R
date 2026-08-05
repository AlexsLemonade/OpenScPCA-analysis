#!/usr/bin/env Rscript
#
# Export SCE reference versions

# sets limit to 48 GB, needed for reading in the qs2 files
mem.maxVSize(48000)

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--input_ref_file"),
    type = "character",
    default = "",
    help = "Path to Seurat version of an OsteoCar object in qs2 format"
  ),
  make_option(
    opt_str = c("--output_sce_file"),
    type = "character",
    help = "Path to RDS file with an SCE version of the reference."
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))

stopifnot("input_ref_file does not exist" = file.exists(opts$input_ref_file))

# load the bigger libraries after passing checks
suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
})

# read input file
osteocar_seurat <- qs2::qs_read(opts$input_ref_file)

# convert to SCE, using `as` to avoid CI error
# the clean up up for memory right away
osteocar_sce <- as.SingleCellExperiment(osteocar_seurat)
rm(osteocar_seurat)
gc()

# remove unneeded items from object
reducedDims(osteocar_sce) <- NULL

readr::write_rds(
  osteocar_sce,
  opts$output_sce_file,
  compress = "gz"
)
