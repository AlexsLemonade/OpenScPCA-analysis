#!/usr/bin/env Rscript


# Load libraries and renv environment ------
library(SingleCellExperiment)
library(optparse)

project_root <- here::here()
renv::load(project_root)


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

  if (is.null(sample_var)) {
    # SCE with single sample
    result_df <- scDblFinder::scDblFinder(
      sce,
      BPPARAM = BiocParallel::MulticoreParam(cores),
      returnType = "table", # return df, not sce
      ...
    )
  } else {
    if (is.null(random_seed)) {
      # See section 1.5.6: https://bioconductor.org/packages/3.19/bioc/vignettes/scDblFinder/inst/doc/scDblFinder.html#usage
      stop("If multiple samples are present in the input SCE object, a random seed must be provided for reproducibility.")
    }

    if (cores == 1) {
      bp <- BiocParallel::SerialParam(RNGseed = random_seed)
    } else {
      bp <- BiocParallel::MulticoreParam(cores, RNGseed = random_seed)
    }
    result_df <- scDblFinder::scDblFinder(
      sce,
      samples = sample_var,
      BPPARAM = bp,
      returnType = "table", # return df, not sce
      ...
    )
  }

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
    default = "",
    help = "The directory containing the input RDS file."
  ),
  make_option(
    "--results_dir",
    type = "character",
    default = "",
    help = "The directory to export TSV file with doublet inferences."
  ),
  make_option(
    "--cores",
    type = "integer",
    default = 4,
    help = "Number of cores to use during scDblFinder inference."
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
if (opts$dataset_name == "") {
  stop("Must provide a dataset name with --dataset_name.")
}
if (!dir.exists(opts$data_dir)) {
  stop("Must provide an input path to data files with --data_dir.")
}
if (opts$results_dir == "") {
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
  run_scdblfinder(cores = opts$cores) |>
  readr::write_tsv(output_tsv_file)
