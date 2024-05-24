#!/usr/bin/env Rscript


# Load libraries and renv environment ------
project_root <- here::here()
renv::load(project_root)

library(SingleCellExperiment)
library(optparse)


# Functions -----------

#' Run scDblFinder on an SCE object
#'
#' @param sce SCE object to run scDblFinder on
#' @param cores Number of cores to specify to BiocParallel::MulticoreParam(
#' @param sample_var Required for SCE object with multiple samples, the SCE colData variable that indicates sample of origin
#' @param random_seed Required for SCE object with multiple samples, a random seed to ensure reproducibility when running on multiple samples
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
      stop("If multiple samples are present in the input SCE object, a random seed must be provided for reproducibility.")
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
    ...
  )

  result_df <- result_df |>
    as.data.frame() |>
    tibble::rownames_to_column("barcodes")

  return(result_df)
}


# Parse options --------

option_list <- list(
  make_option(
    "--dataset_name",
    type = "character",
    help = "Name of dataset to process, where the associated file is expected to be named `{name}_sce.rds`"
  ),
  make_option(
    "--data_dir",
    type = "character",
    help = "The directory containing the input RDS file."
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
    "--random_seed",
    type = "integer",
    default = 2024,
    help = "Random seed. Default is 2024."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Check input arguments
if (is.null(opts$dataset_name)) {
  stop("Must provide a dataset name with --dataset_name.")
}
if (is.null(opts$data_dir)) {
  stop("Must provide an input path for the SCE file with --data_dir.")
}
if (is.null(opts$results_dir)) {
  stop("Must provide an output path for results with --results_dir.")
}

fs::dir_create(opts$results_dir)
set.seed(opts$random_seed)


# Detect doublets and export TSV file with inferences -----
input_sce_file <- file.path(opts$data_dir, glue::glue("{opts$dataset_name}_sce.rds"))
if (!file.exists(input_sce_file)) {
  stop(
    glue::glue("Could not find input file, expected at: `{input_sce_file}`.")
  )
}
output_tsv_file <- file.path(opts$results_dir, glue::glue("{opts$dataset_name}_scdblfinder.tsv"))

readRDS(input_sce_file) |>
  run_scdblfinder(
    cores = opts$cores
  ) |>
  readr::write_tsv(output_tsv_file)
