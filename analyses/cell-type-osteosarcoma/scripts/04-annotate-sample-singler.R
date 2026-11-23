#!/usr/bin/env Rscript
#
# This script performs cell type annotation on one ScPCA SCE sample using a
#  trained SingleR model.
# The script exports an RDS file containing the full SingleR result DataFrame object
#  and a TSV file containing only the annotations and delta.next values.

suppressWarnings({
  suppressPackageStartupMessages({
    library(optparse)
    library(SingleCellExperiment)
  })
})

# Functions --------------------------------------------------------------------

# Annotate an ScPCA SCE object with a trained SingleR model
#
# Inputs:
#   sce: SingleCellExperiment object from ScPCA, with Ensembl ids as row names
#   singler_model: trained SingleR model
#   bp_param: BiocParallel parameter object for SingleR to use
#
# Output: SingleR result DataFrame
annotate_with_singler <- function(sce, singler_model, bp_param = BiocParallel::SerialParam()) {
  # Convert Ensembl ids to gene symbols so they match the model.
  # The OsteoCAR reference uses 10x 2020 (Ensembl v98) symbols, so we map
  #  with the same reference used when training the model.
  # Note that ids that don't map will retain Ensembl ids.
  sce <- rOpenScPCA::sce_to_symbols(sce, reference = "10x2020")
  # Perform the annotation using the provided model and parallelization parameter
  SingleR::classifySingleR(
    sce,
    singler_model,
    BPPARAM = bp_param
  )
}

# Options ----------------------------------------------------------------------

option_list <- list(
  make_option(
    opt_str = c("--sce_file"),
    type = "character",
    default = "",
    help = "Input ScPCA SCE object (processed RDS file) to run SingleR on"
  ),
  make_option(
    opt_str = c("--singler_model_file"),
    type = "character",
    default = "",
    help = "Path to trained SingleR model"
  ),
  make_option(
    opt_str = c("--singler_output_tsv"),
    type = "character",
    default = "",
    help = "Path to output TSV file to save lightweight SingleR annotations"
  ),
  make_option(
    opt_str = c("--singler_output_rds"),
    type = "character",
    default = "",
    help = "Path to output RDS file to save the full SingleR result DataFrame"
  ),
  make_option(
    opt_str = c("--threads"),
    type = "integer",
    default = 4,
    help = "Number of threads for SingleR to use"
  ),
  make_option(
    opt_str = c("--seed"),
    type = "integer",
    default = 2025,
    help = "Random seed"
  )
)

# Parse options and check arguments
opts <- parse_args(OptionParser(option_list = option_list))
stopifnot(
  "sce_file does not exist" = file.exists(opts$sce_file),
  "singler_model_file does not exist" = file.exists(opts$singler_model_file),
  "singler_output_tsv must be provided" = opts$singler_output_tsv != "",
  "singler_output_tsv must end in .tsv or .tsv.gz" =
    stringr::str_detect(opts$singler_output_tsv, "\\.tsv(\\.gz)?$"),
  "singler_output_rds must be provided" = opts$singler_output_rds != "",
  "singler_output_rds must end in .rds" =
    stringr::str_detect(opts$singler_output_rds, "\\.rds$")
)
set.seed(opts$seed)

# Set up parallelization parameters for SingleR
if (opts$threads == 1) {
  bp_param <- BiocParallel::SerialParam()
} else {
  bp_param <- BiocParallel::MulticoreParam(opts$threads)
}

# Make sure output directories exist
fs::dir_create(dirname(opts$singler_output_tsv))
fs::dir_create(dirname(opts$singler_output_rds))

# Annotate ---------------------------------------------------------------------
#  Read in the SingleR model and the SCE object
singler_model <- readRDS(opts$singler_model_file)
sce <- readRDS(opts$sce_file)
#  Annotate the SCE object with the SingleR model
singler_result <- annotate_with_singler(sce, singler_model, bp_param)

# Export -----------------------------------------------------------------------
# Extract TSV to save separately
singler_df <- singler_result |>
  # keep only labels, delta.next, and pruned.labels
  purrr::discard_at("scores") |>
  as.data.frame() |>
  tibble::rownames_to_column("barcodes")

# Export the results to RDS and TSV files
readr::write_rds(singler_result, opts$singler_output_rds)
readr::write_tsv(singler_df, opts$singler_output_tsv)
