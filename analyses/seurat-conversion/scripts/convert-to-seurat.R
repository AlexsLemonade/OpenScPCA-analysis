#!/usr/bin/env Rscript

# This script takes a directory containing ScPCA-formatted SingleCellExperiment
# objects (filtered or processed only) as `.rds` files and converts each object
# to a Seurat object. The Seurat object is then saved as an `rds` file in the
# specified output directory.


library(optparse)

# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-i", "--input_dir"),
    type = "character",
    help = "Path where SCE rds files are located. Will search recursively for files ending in `_processed.rds`, `_filtered.rds`, or `_unfiltered.rds.`",
  ),
  make_option(
    opt_str = c("-o", "--output_dir"),
    type = "character",
    help = "Path where Seurat rds files will be saved.",
  ),
  make_option(
    opt_str = c("--use_ensembl"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Use gene symbols for Seurat object. Default is to use Ensembl IDs."
  ),
  make_option(
    opt_str = c("--dedup_method"),
    type = "character",
    default = "sum",
    help = "Method to use for deduplication. Default is 'sum', with the alternative 'unique'."
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

stopifnot(
  "The input directory does not exist" = file.exists(opts$input_dir),
  "--dedup_method must be one of 'sum' or 'unique'" = opts$dedup_method %in% c("sum", "unique")
)

# Load libraries ---------------------------------------------------------------

# get input file list with relative paths
sce_paths <- list.files(
  opts$input_dir,
  pattern = "_(processed|filtered).rds$",
  recursive = TRUE
)

if (length(sce_paths) == 0) {
  warning("No SCE files found in input directory: ", opts$input_dir)
}

# There are some annoying warning when the SCE package loads (even implicitly)
# let's suppress them
suppressPackageStartupMessages({
  suppressWarnings(
    library(SingleCellExperiment)
  )
})

# convert each file
sce_paths |> purrr::walk(\(path){
  message("Converting ", path, " to Seurat format")
  sce_file <- file.path(opts$input_dir, path)

  # generate output file name and create path if it doesn't exist
  seurat_file <- file.path(
    opts$output_dir,
    stringr::str_replace(path, ".rds$", "_seurat.rds")
  )
  fs::dir_create(dirname(seurat_file))

  sce <- readRDS(sce_file)
  seurat_obj <- rOpenScPCA::sce_to_seurat(
    sce,
    use_symbols = !opts$use_ensembl,
    dedup_method = opts$dedup_method
  )
  saveRDS(seurat_obj, seurat_file)

})
