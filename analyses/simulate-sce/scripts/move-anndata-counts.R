#!usr/bin/env Rscript

# Script to move the logcounts layer to the X slot in an anndata object

library(basilisk)
library(optparse)

# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-d", "--dir"),
    type = "character",
    help = "Path where AnnData files are located. Will search recursively for files ending in `_processed_*.h5ad",
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# get processed files

anndata_files <- list.files(
  opts$dir,
  pattern = "_processed_.+\\.h5ad$",
  recursive = TRUE,
  full.names = TRUE
)


# load basilisk
proc <- basiliskStart(env = zellkonverter::zellkonverterAnnDataEnv(), testload = 'anndata')
on.exit(basiliskStop(proc))


# run a function in the basilisk environment to move elements of files
basiliskRun(proc, fun = function(files){
  adata <- reticulate::import("anndata")
  for(afile in files){
    h5ad <- adata$read_h5ad(afile)
    if (!is.null(h5ad$layers$get("logcounts"))){
      h5ad$raw <- h5ad
      h5ad$X <- h5ad$layers["logcounts"]
      h5ad$uns$X_name <- "logcounts"
      h5ad$layers["logcounts"] <- NULL
      h5ad$write_h5ad(afile, compression = "gzip")
    }
  }
}, files = anndata_files)
