#!/usr/bin/env Rscript
#
# This is a helper script to extract all cell ids from the full v2 NBAtlas object `seuratObj_NBAtlas_share_v20241203.rds` from
# this Mendeley record: https://data.mendeley.com/datasets/yhcf6787yp/3
# Any cell ids _not_ present in this object will not be included in our converted references, created with `00b_convert-nbatlas.R`


library(optparse)
library(Seurat)

option_list <- list(
  make_option(
    opt_str = c("--nbatlas_file"),
    type = "character",
    default = "scratch/seuratObj_NBAtlas_share_v20241203.rds",
    help = "Path to Seurat version of an NBAtlas object"
  ),
  make_option(
    opt_str = c("--cell_id_file"),
    type = "character",
    default = "scratch/nbtatlas-ids.txt.gz",
    help = "Path to output file to store cell ids"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))
stopifnot("nbatlas_file does not exist" = file.exists(opts$nbatlas_file))

# Extract ids and save
readRDS(opts$nbatlas_file) |>
  colnames() |>
  readr::write_lines(opts$cell_id_file)
