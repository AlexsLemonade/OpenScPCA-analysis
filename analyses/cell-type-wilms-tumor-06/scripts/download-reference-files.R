#!/usr/bin/env Rscript

# This script downloads files used as part of label transfer:
# - The fetal kidney atlas (Stewart) RDS file from CELLxGENE. This is saved in <output_dir>/fetal_kidney.rds
# - The homologs.rds file from Azimuth. This is saved in <output_dir>/homologs.rds
#
# USAGE:
# Rscript download-fetal-reference.R \
#   --output_dir scratch

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--kidney_url"),
    type = "character",
    default = "https://datasets.cellxgene.cziscience.com/40ebb8e4-1a25-4a33-b8ff-02d1156e4e9b.rds",
    help = "The URL of the fetal kidney atlas from CELLxGENE"
  ),
  make_option(
    opt_str = c("--homologs_url"),
    type = "character",
    default = "https://seurat.nygenome.org/azimuth/references/homologs.rds",
    help = "The URL of the homologs.rds file from Azimuth. It will be saved as homologs.rds."
  ),
  make_option(
    opt_str = c("-d", "--output_dir"),
    type = "character",
    default = "scratch",
    help = "Output directory for the fetal kidney atlas rds file, relative to the current directory. It will be saved as fetal_kidney.rds."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))


dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

kidney_file <- file.path(
  opts$output_dir,
  "fetal_kidney.rds"
)
homologs_file <- file.path(
  opts$output_dir,
  "homologs.rds"
)


if (!file.exists(kidney_file)) {
  download.file(url = opts$kidney_url, destfile = kidney_file)
}

if (!file.exists(homologs_file)) {
  download.file(url = opts$homologs_url, destfile = homologs_file)
}
