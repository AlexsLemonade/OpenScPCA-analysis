#!/usr/bin/env Rscript

# Download the fetal kidney dataset and create a reference for use with
# Azimuth
#
# USAGE:
# Rscript download-build-kidney-reference.R \
#   --url https://datasets.cellxgene.cziscience.com/40ebb8e4-1a25-4a33-b8ff-02d1156e4e9b.rds \
#   --output_dir ../results/references \
#   --seed 2024
#

library(optparse)
library(Seurat)
library(Azimuth)

# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-u", "--url"),
    type = "character",
    default = "https://datasets.cellxgene.cziscience.com/40ebb8e4-1a25-4a33-b8ff-02d1156e4e9b.rds",
    help = "The URL of the fetal kidney atlas from CELLxGENE"
  ),
  make_option(
    opt_str = c("-d", "--output_dir"),
    type = "character",
    default = "results/references",
    help = "Output directory for the Azimuth reference, relative to your current directory"
  ),
  make_option(
    opt_str = c("-s", "--seed"),
    type = "integer",
    default = 12345,
    help = "Seed passed to set.seed()"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# Download data ----------------------------------------------------------------

project_root <- rprojroot::find_root(rprojroot::is_git_root)
path_to_data <- file.path(
  project_root,
  "analyses",
  "cell-type-wilms-tumor-06",
  "scratch",
  "fetal_kidney.rds"
)
download.file(url = opts$url, destfile = path_to_data)

# Read in data -----------------------------------------------------------------

seurat <- readRDS(path_to_data)

# Transform and dimension reduction --------------------------------------------

set.seed(opts$seed)

options(future.globals.maxSize = 891289600000)
s <- SCTransform(
  seurat,
  verbose = FALSE,
  method = "glmGamPoi",
  conserve.memory = TRUE
)
s <- RunPCA(s, npcs = 50, verbose = FALSE)
s <- RunUMAP(s, dims = 1:50, verbose = FALSE, return.model = TRUE)

# Create reference -------------------------------------------------------------

options(future.globals.maxSize = 891289600000)
fetal_kidney <- AzimuthReference(
  s,
  refUMAP = "umap",
  refDR = "pca",
  refAssay = "SCT",
  dims = 1:50,
  k.param = 31,
  plotref = "umap",
  plot.metadata = NULL,
  ori.index = NULL,
  colormap = NULL,
  assays = NULL,
  metadata = c("compartment", "cell_type"),
  reference.version = "0.0.0",
  verbose = FALSE
)

# Save reference ---------------------------------------------------------------

# Create directory if it doesn't exist yet
dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)

# Save annoy index
SaveAnnoyIndex(
  object = fetal_kidney[["refdr.annoy.neighbors"]],
  file = file.path(opts$output_dir, "idx.annoy")
)
# Save reference
saveRDS(object = fetal_kidney, file = file.path(opts$output_dir, "ref.Rds"))