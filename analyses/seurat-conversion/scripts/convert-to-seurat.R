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
    help = "Path where SCE rds files are located. Will search recursively for files ending in `_processed.rds` or  `_filtered.rds`.",
  ),
  make_option(
    opt_str = c("-o", "--output_dir"),
    type = "character",
    help = "Path where Seurat rds files will be saved.",
  ),
  make_option(
    opt_str = c("--use_ensembl"),
    dest = "use_symbols",
    type = "logical",
    default = TRUE,
    action = "store_false",
    help = "Use Ensembl ids for Seurat object. Default is to use gene symbols."
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

# get input file list with relative paths
sce_paths <- list.files(
  opts$input_dir,
  pattern = "_(processed|filtered).rds$",
  recursive = TRUE
)

if (length(sce_paths) == 0) {
  warning("No SCE files found in input directory: ", opts$input_dir)
}

# There are several warnings when the SCE package loads (even implicitly)
# Let's suppress them
suppressPackageStartupMessages({
  suppressWarnings(
    library(SingleCellExperiment)
  )
})

# convert each file
sce_paths |> purrr::walk(\(path){
  message("Converting ", path, " to Seurat format")
  sce_file <- file.path(opts$input_dir, path)

  # generate output file name
  seurat_file <- file.path(
    opts$output_dir,
    stringr::str_replace(path, ".rds$", "_seurat.rds")
  )
  # create output subdirectory if it doesn't exist
  fs::dir_create(dirname(seurat_file))

  # perform processing
  sce <- readRDS(sce_file)
  seurat_obj <- rOpenScPCA::sce_to_seurat(
    sce,
    use_symbols = opts$use_symbols,
    dedup_method = opts$dedup_method
  )
  saveRDS(seurat_obj, seurat_file)

})
