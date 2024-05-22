#!/usr/bin/env Rscript

library(SingleCellExperiment)

# Setup ------

# Define directories
module_base <- file.path(
  rprojroot::find_root(rprojroot::is_git_root),
  "analyses",
  "doublet-detection"
)

data_dir <- file.path(module_base, "scratch", "benchmark_datasets")
result_dir <- file.path(module_base, "results", "benchmark_results")
fs::dir_create(result_dir)


# Functions -----------

#' Run scDblFinder on an SCE object
#'
#' @param sce SCE object to run scDblFinder on
#' @param cores Number of cores to specify to BiocParallel::MulticoreParam(
#' @param sample_var Required for multiplexed data, the SCE colData variable that indicates sample of origin
#' @param random_seed Required for multiplexed data, a random seed to ensure reproducibility when running on multiplexed data.
#' @param ... Additional arguments to pass to scDblFinder::scDblFinder
#'
#' @return SCE object with scDblFinder inferences
run_scdblfinder <- function(sce,
                            cores = 4,
                            sample_var = NULL,
                            random_seed = NULL,
                            ...) {

  if (is.null(sample_var)) {
    # not multiplexed
    sce <- scDblFinder::scDblFinder(sce,
                                    BPPARAM = BiocParallel::MulticoreParam(cores),
                                    ...)
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
    sce <- scDblFinder::scDblFinder(sce,
                                    samples = sample_var,
                                    BPPARAM = bp,
                                    ...)
  }

  return(sce)
}

# Infer doublets --------------

# Benchmark dataset names
datanames <- c("hm-6k", "pbmc-1B-dm", "pdx-MULTI", "HMEC-orig-MULTI")

# Detect doublets and export SCE file with inferences
datanames |>
  purrr::walk(
    \(dataname) {
      input_sce_file <- file.path(data_dir, glue::glue("{dataname}_sce.rds"))
      output_sce_file <- file.path(result_dir, glue::glue("{dataname}_sce.rds"))

      readRDS(input_sce_file) |>
        run_scdblfinder() |>
        readr::write_rds(output_sce_file)
    }
  )
