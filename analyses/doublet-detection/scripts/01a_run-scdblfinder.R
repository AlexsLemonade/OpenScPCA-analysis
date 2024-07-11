#!/usr/bin/env Rscript

# This script runs scDblFinder on an SCE file and exports a TSV file of results.


# Load libraries and renv environment ------
project_root <- rprojroot::find_root(rprojroot::is_renv_project)
renv::load(project_root)

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(optparse)
})

# Functions -----------

#' Run scDblFinder on an SCE object
#'
#' @param sce SCE object to run scDblFinder on
#' @param cores Number of cores to specify to BiocParallel::MulticoreParam()
#' @param sample_var Required for SCE object with multiple samples, the SCE colData variable that indicates sample of origin.
#'   This option should _not_ be used for a multiplexed library which has not been demultiplexed, but is suitable for merged objects.
#' @param random_seed Required for SCE object with multiple samples, a random seed to ensure reproducibility when running on multiple samples.
#' @param ... Additional arguments to pass to scDblFinder::scDblFinder
#'
#' @return Data frame object with scDblFinder inferences
run_scdblfinder <- function(sce,
                            cores = 4,
                            sample_var = NULL,
                            random_seed = NULL,
                            ...) {

  # first check multiple samples, which requires a random seed
  if (!is.null(sample_var) && is.null(random_seed)) {
    if (!sample_var %in% colnames(colData(sce))) {
      stop(
        glue::glue("The provided sample variable {sample_var} is not present in the input SCE's colData slot.")
      )
    }
    if (is.null(random_seed)) {
      # See section 1.5.6: https://bioconductor.org/packages/3.19/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html#usage
      stop("If multiple samples are present in the input SCE object (e.g., the object contains multiple merged libraries), a random seed must be provided for reproducibility.")
    }
  }

  # set up cores
  if (cores == 1) {
    bp <- BiocParallel::SerialParam(RNGseed = random_seed)
  } else {
    bp <- BiocParallel::MulticoreParam(cores, RNGseed = random_seed)
  }

  # Run doublet finder
  result_df <- scDblFinder::scDblFinder(
    sce,
    samples = sample_var, # Default is NULL so use whatever was provided by the user or default
    BPPARAM = bp,
    returnType = "table", # return df, not sce
    verbose = FALSE,
    ...
  )

  result_df <- result_df |>
    as.data.frame() |>
    tibble::rownames_to_column("barcodes") |>
    # remove artifical doublets
    dplyr::filter(!stringr::str_starts(barcodes, "rDbl"))

  return(result_df)
}


# Parse options --------

option_list <- list(
  make_option(
    "--input_sce_file",
    type = "character",
    help = "Path to input SCE file in RDS format."
  ),
  make_option(
    "--results_dir",
    type = "character",
    help = "The directory to export TSV file with doublet inferences."
  ),
  make_option(
    "--cores",
    type = "integer",
    default = 4,
    help = "Number of cores to use during scDblFinder inference. Only used when there are multiple samples in the SCE."
  ),
  make_option(
    "--sample_var",
    type = "character",
    help = "If multiple samples are present in the SCE, the colData column name with sample ids. This option should be used for merged objects."
  ),
  make_option(
    "--random_seed",
    type = "integer",
    default = 2024,
    help = "Random seed. Default is 2024."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check input arguments and set up files, directories ------
if (!file.exists(opts$input_sce_file)) {
    glue::glue("Could not find input SCE file, expected at: `{opts$input_sce_file}`.")
}
if (is.null(opts$results_dir)) {
  stop("Must provide an output path for results with --results_dir.")
}

output_tsv_file <- file.path(
  opts$results_dir,
  stringr::str_replace(
    basename(opts$input_sce_file),
    ".rds",
    "_scdblfinder.tsv"
  )
)

fs::dir_create(opts$results_dir)
set.seed(opts$random_seed)

# Detect doublets and export TSV file with inferences -----

readRDS(opts$input_sce_file) |>
  run_scdblfinder(
    cores = opts$cores,
    # only used if there are multiple samples, e.g. if this is processing a merged object
    sample_var = opts$sample_var, # will be NULL if not provided as input argument
    random_seed = opts$random_seed
  ) |>
  readr::write_tsv(output_tsv_file)
