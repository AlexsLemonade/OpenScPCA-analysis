#!/usr/bin/env Rscript

# This script runs four doublet detection approaches on a set of selected samples, for non-multiplexed projects
# Methods include scDblFinder, scds::cxds, scds::bcds, and scds::hybrid.
# Seven libraries (where possible) of varying library sizes per project are identified, and doublets are detected in each.
# The script outputs two TSVs to `../results`:
#  - TSV of runtimes in seconds for each library
#  - TSV of doublet detection results for each library
# SCEs with doublet inferences are also saved in `../scratch/benchmark-sces`


# Load libraries -------------------
library(scds)
library(scDblFinder)
library(optparse)

# Define functions to detect doublets using four methods ---------

capture_time <- function(start_time, end_time) {
  # Helper function to determine runtime given `start_time` and `end_time` values as calculated with `Sys.time()`
  # This function returns a numeric value, representing runtime in seconds
  (end_time - start_time) |>
    as.numeric(units = "secs") |>
    round(2)
}


benchmark_doublets <- function(input_sce_file,
                               library_id,
                               scratch_dir) {
  # Run doublet detection methods on a given input SCE file and its library_id
  # This function exports two TSVs into the scratch directory containing runtime and doublet results, to support faster script re-execution
  # This function returns a list of two data frames:
  #  - `results` contains the results from each doublet detection methods
  #  - `runtime` contains the runtime, in seconds, for each doublet detection method

  if (!(file.exists(input_sce_file))) {
    stop(glue::glue("Missing SCE file for {library_id}"))
  }

  library_results_tsv <- file.path(
    scratch_dir,
    glue::glue("{library_id}_benchmark-results.tsv")
  )
  library_runtime_tsv <- file.path(
    scratch_dir,
    glue::glue("{library_id}_benchmark-runtime.tsv")
  )

  if (!file.exists(library_results_tsv) & !file.exists(library_runtime_tsv)) {

    sce <- readRDS(input_sce_file)

    # scDblFinder, with a bonus cxds variant in the `cxds_score` column
    # https://github.com/plger/scDblFinder/blob/99b947a49a4f3f3bd1a2203fe2296a6dc143fe11/R/misc.R#L266-L274
    #
    # Note further, this warns "scDblFinder might not work well with very low numbers of cells." when <100 cells
    start.time <- Sys.time()
    sce <- scDblFinder(
      sce,
      clusters = sce$clusters,
      BPPARAM = BiocParallel::MulticoreParam(opts$cores)
    )
    end.time <- Sys.time()
    metadata(sce)$scdbl_time <- capture_time(start.time, end.time)
    scdbl_df <- tibble::tribble(
      ~barcode,
      sce$barcodes,
    )

    # hybrid
    start.time <- Sys.time()
    sce <- cxds_bcds_hybrid(sce)
    end.time <- Sys.time()
    metadata(sce)$hybrid_time <- capture_time(start.time, end.time)

    # cxds
    start.time <- Sys.time()
    sce <- cxds(sce)
    end.time <- Sys.time()
    metadata(sce)$cxds_time <- capture_time(start.time, end.time)

    # bcds
    start.time <- Sys.time()
    sce <- bcds(sce)
    end.time <- Sys.time()
    metadata(sce)$bcds_time <- capture_time(start.time, end.time)

    # Save the SCE with doublet results
    readr::write_rds(sce, file.path(
      scratch_dir,
      glue::glue("{library_id}-doublet-benchmark.rds")
    ))


  results_df <- tibble::tribble(
    ~scpca_library_id, ~barcode,     ~scdbl_class,          ~scdbl_cxds,               ~sbdbl_score,           ~cxds_score,    ~bcds_score,    ~hybrid_score,
    library_id,       sce$barcodes, sce$scDblFinder.class, sce$scDblFinder.cxds_score, sce$scDblFinder.score, sce$cxds_score, sce$bcds_score, sce$hybrid_score
  ) |>
    tidyr::unnest(cols = -scpca_library_id)

  runtimes_df <- tibble::tribble(
    ~scpca_library_id, ~scdbl_s,                  ~cxds_s,                 ~bcds_s,                 ~hybrid_s,
    library_id,        metadata(sce)$scdbl_time, metadata(sce)$cxds_time, metadata(sce)$bcds_time, metadata(sce)$hybrid_time
  )

  # memory purge
  rm(sce)

  # save them into scratch
  readr::write_tsv(results_df, library_results_tsv)
  readr::write_tsv(runtimes_df, library_runtime_tsv)

  } else {

    # If sample was already processed, just read its results in
    results_df <- readr::read_tsv(library_results_tsv)
    runtimes_df <- readr::read_tsv(library_runtime_tsv)
  }

  return(
    list(
    "results" = results_df,
    "runtimes" = runtimes_df
  ))

}

# Parse arguments ----------------------------------------

option_list <- list(
  make_option(
    c("-c", "--cores"),
    type = "integer",
    default = 4,
    help = "Number of cores to use for parallel processing"
  )
)
opts <- parse_args(OptionParser(option_list = option_list))


# Define paths and set seed ------------------
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
data_dir <- file.path(repository_base, "data", "current")
module_base <- file.path(repository_base, "analyses", "doublet-detection")

# directory to save doublet detection results
output_dir <- file.path(module_base, "results", "benchmark-results")
fs::dir_create(output_dir)

# directory to save temporary files
scratch_dir <- file.path(module_base, "scratch", "benchmark-sces")
fs::dir_create(scratch_dir)

# output files
benchmark_runtime_tsv <- file.path(output_dir, "benchmark_runtimes.tsv")
benchmark_results_tsv <- file.path(output_dir, "benchmark_results.tsv")

set.seed(11)

# Read in all sample metadata -----------------
metadata_files <- list.files(
  data_dir,
  pattern = "single_cell_metadata.tsv",
  full.names = TRUE,
  recursive = TRUE
)
names(metadata_files) <- stringr::str_extract(metadata_files, "SCPCP\\d{6}")

# Exclude multiplexed, for now
metadata_files <- metadata_files[names(metadata_files) != "SCPCP000009"]

# data frame of all library diagnoses and processed cell counts
sample_df <- metadata_files |>
  purrr::imap(
    \(file, project_id) {
      readr::read_tsv(file) |>
        dplyr::select(scpca_sample_id,
                      scpca_library_id,
                      diagnosis,
                      processed_cells) |>
        dplyr::mutate(project_id = project_id) |>
        dplyr::select(project_id, everything())

    }
  ) |>
  dplyr::bind_rows()

# Identify libraries to use ------------------------

# We'll keep 7 libraries per project (or all libraries, if < 10 total)
# across a range of library sizes - rank each and select bottom 2, middle 3, top 2

min_cells <- 50

libraries_to_benchmark <- sample_df |>
  # keep only libraries with >= min_cells cells
  dplyr::filter(processed_cells >= min_cells) |>
  dplyr::group_by(project_id) |>
  dplyr::mutate(n_samples = dplyr::n()) |>
  dplyr::arrange(processed_cells) |>
  dplyr::slice(c(
    1:2, # bottom 2
    (floor(n_samples/2) - 1):(floor(n_samples/2) + 1), # roughly middle 3
    (n_samples - 1):n_samples # top 2
  )) |>
  dplyr::ungroup() |>
  # remove duplicates for projects with <10 libraries
  dplyr::distinct() |>
  dplyr::select(-n_samples) |>
  # add column with file path
  dplyr::mutate(file_path = file.path(
    data_dir, project_id, scpca_sample_id, glue::glue("{scpca_library_id}_processed.rds")
  ))

# Print samples as a comma-sep list for data download
#print(paste0(libraries_to_benchmark$scpca_sample_id, collapse = ","))

# Run doublet detection methods and export results ------------------------------------

sce_files <- libraries_to_benchmark$file_path
names(sce_files) <- libraries_to_benchmark$scpca_library_id

doublet_output <- sce_files|>
  purrr::imap(benchmark_doublets, scratch_dir)

runtime_df <- doublet_output |>
  purrr::map(\(x) x$runtimes) |>
  dplyr::bind_rows()

result_df <- doublet_output |>
  purrr::map(\(x) x$results) |>
  dplyr::bind_rows()


# Export runtimes, combined with some sample metadata
libraries_to_benchmark |>
  dplyr::select(-file_path) |>
  dplyr::inner_join(
    runtime_df,
    by = "scpca_library_id"
  ) |>
  readr::write_tsv(benchmark_runtime_tsv)

# Export results
readr::write_tsv(
  result_df,
  benchmark_results_tsv
)


