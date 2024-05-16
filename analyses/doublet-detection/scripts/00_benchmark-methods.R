# This script runs four doublet detection approaches on a set of selected samples, for non-multiplexed projects
# Seven libraries (where possible) of varying library sizes per project are identified, and doublets are detected in each.
# Outputs a TSV of run times in seconds for each library and SCE objects with doublet inferences in the colData for subsequent analysis.


# Load libraries -------------------
library(scds)
library(scDblFinder)
library(optparse)

# Define function to detect doublets using four methods ---------

benchmark_doublets <- function(input_sce_file, sample_id, output_dir) {

  if (!(file.exists(input_sce_file))) {
    stop(glue::glue("Missing SCE file for {sample_id}"))
  }

  output_sce_path <- file.path(
    output_dir,
    glue::glue("{sample_id}-doublet-benchmark.rds")
  )
  if (file.exists(output_sce_path)) {
    sce <- readRDS(output_sce_path)
  } else {

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
    metadata(sce)$scdbl_time <- (end.time - start.time) |>
      as.numeric(units = "secs") |>
      round(2)

    # hybrid
    start.time <- Sys.time()
    sce <- cxds_bcds_hybrid(sce)
    end.time <- Sys.time()
    metadata(sce)$hybrid_time <- (end.time - start.time) |>
      as.numeric(units = "secs") |>
      round(2)

    # csds
    start.time <- Sys.time()
    sce <- cxds(sce)
    end.time <- Sys.time()
    metadata(sce)$cxds_time <- (end.time - start.time) |>
      as.numeric(units = "secs") |>
      round(2)

    # bsds
    start.time <- Sys.time()
    sce <- bcds(sce)
    end.time <- Sys.time()
    metadata(sce)$bsds_time <- (end.time - start.time) |>
      as.numeric(units = "secs") |>
      round(2)

    # Save the SCE with doublet results
    readr::write_rds(sce, file.path(
      output_dir,
      glue::glue("{sample_id}-doublet-benchmark.rds")
    ))
  }

  # Create and return tibble of runtimes
  return(
    tibble::tribble(
      ~scpca_sample_id, ~scdblfinder_s, ~cxds_s,       ~bsds_s,       ~hybrid_s,
      sample_id,        metadata(sce)$scdbl_time, metadata(sce)$cxds_time, metadata(sce)$bsds_time, metadata(sce)$hybrid_time
    )
  )

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


# Define paths ------------------
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
data_dir <- file.path(repository_base, "data", "current")
module_base <- file.path(repository_base, "analyses", "doublet-detection")

# directory to save SCEs with doublet detection results
output_dir <- file.path(module_base, "scratch", "benchmark-sces")
fs::dir_create(output_dir)

# output TSV file with benchmarking times
benchmark_tsv <- file.path(module_base,
                           "results",
                           glue::glue("benchmark-times_{opts$cores}-cores.tsv"))

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
                      #unfiltered_cells,
                      #filtered_cells = filtered_cell_count,
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
    (floor(n_samples/2) - 1):(floor(n_samples/2) + 1), # roughly middle 3,
    (n_samples -2):n_samples # top 2
  )
) |>
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

# Run initial doublet detection methods ------------------------------------

sample_files <- libraries_to_benchmark$file_path
names(sample_files) <- libraries_to_benchmark$scpca_sample_id

times_df <- sample_files |>
  purrr::imap(benchmark_doublets, output_dir) |>
  # combine runtime results with sample metadata
  dplyr::bind_rows()

libraries_to_benchmark |>
  dplyr::select(-file_path) |>
  dplyr::left_join(
    times_df,
    by = "scpca_sample_id"
  ) |>
  readr::write_tsv(benchmark_tsv)
