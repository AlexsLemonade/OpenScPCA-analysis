#!/usr/bin/env Rscript

# This script prepares references uses for label transfer, which will be performed with
# code adapted from Azimuth. Therefore, references should be in a format expected by Azimuth, and
# then they need to be further prepared for use with our adapted code.
#
# Two references are prepared:
# - A fetal kidney atlas (Stewart) obtained from CELLxGENE. This is saved to <output_dir>/stewart_formatted_ref.rds
# - A fetal organ reference (Cao) obtained from Azimuth. This is saved to <output_dir>/cao_formatted_ref.rds
#
# USAGE:
# Rscript prepare-fetal-references.R \
#   --kidney_ref_file ../scratch/fetal_kidney.rds \
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
    opt_str = c("--kidney_ref_file"),
    type = "character",
    default = "scratch/fetal_kidney.rds",
    help = "The relative path from the current directory to the fetal kidney atlas downloaded from CELLxGENE"
  ),
  make_option(
    opt_str = c("-d", "--output_dir"),
    type = "character",
    default = "results/references",
    help = "Output directory for the Azimuth references, relative to the current directory"
  ),
  make_option(
    opt_str = c("-s", "--seed"),
    type = "integer",
    default = 12345,
    help = "Seed passed to set.seed()"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# Functions ----------------------------------------------

#' Prepare an Azimuth reference for label transfer
#'
#' This function prepares an Azimuth-formatted reference for label transfer with our
#' code which is adapted from Azimuth::RunAzimuth
#'
#' @param reference The Azimuth reference object to prepare
#' @param annotation_levels A vector of annotation columns present in the Azimuth object
#'
#' @return A list of objects needed for label transfer, including:
#' - reference
#' - reference_dims
#' - refdata
#' - annotation_levels
prepare_azimuth_reference <- function(reference, annotation_levels) {
  # Update column names from `refdr_<1-100>` --> `refDR_<1-100>`
  # Source: https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R#L72-L80
  key.pattern <- "^[^_]*_"
  new.colnames <- gsub(
    pattern = key.pattern,
    replacement = Key(reference[["refDR"]]),
    x = colnames(
      Loadings(
        object = reference[["refDR"]],
        projected = FALSE
      )
    )
  )
  colnames(Loadings(
    object = reference[["refDR"]],
    projected = FALSE
  )) <- new.colnames

  # Determine number of dimensions
  # Source: https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R#L81
  dims <- as.double(slot(reference, "neighbors")$refdr.annoy.neighbors@alg.info$ndim)

  # Create the corresponding refdata object, which uses the reference's annotation levels
  # Source: https://github.com/satijalab/azimuth/blob/243ee5db80fcbffa3452c944254a325a3da2ef9e/R/azimuth.R#L152-L155
  reference_rownames <- rownames(reference)
  refdata <- lapply(X = annotation_levels, function(x) {
    reference[[x, drop = TRUE]]
  })
  names(refdata) <- annotation_levels

  # Create and return list of reference components
  reference_list <- list(
    reference = reference,
    dims = dims,
    annotation_levels = annotation_levels,
    refdata = refdata
  )

  return(reference_list)
}

# Define paths -------------------------------------------


# Create output directory and define paths for final reference files
dir.create(opts$output_dir, recursive = TRUE, showWarnings = FALSE)
stewart_ref_file <- file.path(
  opts$output_dir,
  "stewart_formatted_ref.rds"
)
cao_ref_file <- file.path(
  opts$output_dir,
  "cao_formatted_ref.rds"
)

# Define the annotation levels to use for each reference ----------------
stewart_annotation_levels <- c("compartment", "cell_type")
cao_annotation_levels <- c("annotation.l1", "annotation.l2", "organ")


# Prepare Stewart (fetal kidney) reference ------------------------------

# Read in data
if (!file.exists(opts$kidney_ref_file)) {
  stop("The kidney reference file does not exist. Make sure to run download-reference-files.R to download it first.")
}
seurat <- readRDS(opts$kidney_ref_file)

# Transform and dimension reduction
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

# Create reference
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
  metadata = stewart_annotation_levels,
  reference.version = "0.0.0",
  verbose = FALSE
)

# format for label transfer
stewart_ref_list <- prepare_azimuth_reference(fetal_kidney, stewart_annotation_levels)

# export formatted reference
saveRDS(stewart_ref_list, stewart_ref_file)


# Prepare Cao (full fetal organ) reference ------------------------------

# Load in the reference, keeping only the $map portion
fetus_ref <- SeuratData::LoadData("fetusref", type = "azimuth")$map


# format for label transfer
cao_ref_list <- prepare_azimuth_reference(fetus_ref, cao_annotation_levels)

# export formatted reference
saveRDS(cao_ref_list, cao_ref_file)
