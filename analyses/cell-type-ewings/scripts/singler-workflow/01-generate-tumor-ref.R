#!/usr/bin/env Rscript

# This script is used to create a reference SCE object that contains all high-confident
# tumor cells for all samples in SCPCP000015
# High confident tumor cells are those that are classified as "Tumor" in the `aucell-annotation.sh` workflow

project_root <- here::here()

library(optparse)
library(SingleCellExperiment)

option_list <- list(
  make_option(
    opt_str = c("--data_dir"),
    type = "character",
    default = file.path("..", "..", "data", "current", "SCPCP000015"),
    help = "Path to folder containing all processed SCE files."
  ),
  make_option(
    opt_str = c("--auc_results_dir"),
    type = "character",
    default = file.path(project_root, "results", "aucell_annotation"),
    help = "Path to folder containing output from `aucell-annotation.sh`."
  ),
  make_option(
    opt_str = c("--output_reference_file"),
    type = "character",
    default = file.path(project_root, "scratch", "tumor-ref-singler.rds"),
    help = "Path to file where RDS file will be saved containing all tumor cells."
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# Define functions -------------------------------------------------------------

# subset sce to only have tumor cells present
# if no tumor cells, then don't return anything
subset_tumor_cells <- function(sce_file,
                               annotation_file) {
  # pull out tumor cells from annotation file
  tumor_cells <- readr::read_tsv(annotation_file) |>
    dplyr::filter(auc_classification == "Tumor") |>
    dplyr::pull(barcodes)

  if (length(tumor_cells) == 0) {
    # exit if no tumor cells
    return(NULL)
  }

  # read in sce
  sce <- readr::read_rds(sce_file)

  # pull out library id and participant id
  library_id <- metadata(sce)$library_id
  participant_id <- metadata(sce)$sample_metadata |>
    dplyr::pull(participant_id)

  # filter to only contain tumor cells
  filtered_sce <- sce[, tumor_cells]

  # add tumor labels to coldata
  filtered_sce$library_id <- library_id
  filtered_sce$participant_id <- participant_id
  filtered_sce$ref_tumor_label <- glue::glue("tumor-{library_id}")

  return(filtered_sce)
}

# Set up -----------------------------------------------------------------------

# Define file paths to annotation files
annotation_file_ext <- "_auc-classifications.tsv"
annotation_files <- list.files(opt$auc_results_dir, pattern = annotation_file_ext, recursive = TRUE, full.names = TRUE)

# name file list with library ids
annotation_library_ids <- basename(annotation_files) |>
  stringr::str_remove(annotation_file_ext)
names(annotation_files) <- annotation_library_ids


# Define file paths to sce files
sce_file_ext <- "_processed.rds"
sce_files <- list.files(opt$data_dir, pattern = sce_file_ext, recursive = TRUE, full.names = TRUE)

# name file list with library ids
sce_library_ids <- basename(sce_files) |>
  stringr::str_remove(sce_file_ext)
names(sce_files) <- sce_library_ids

# make sure annotation files are in the same order as sce files
annotation_files <- annotation_files[order(annotation_files, sce_library_ids)]

# Create merged tumor cell ref -------------------------------------------------

# list of tumor cell sce
sce_list <- purrr::map2(sce_files, annotation_files, subset_tumor_cells)

# remove any objects that don't have tumor cells
sce_absent <- sce_list |>
  purrr::map_lgl(is.null)

remaining_sce_list <- sce_list[which(!sce_absent)]

# merge sces into a single reference
merged_ref_sce <- scpcaTools::merge_sce_list(
  sce_list = remaining_sce_list,
  cell_id_column = "barcodes",
  retain_coldata_cols = c("library_id", "participant_id", "ref_tumor_label")
)


# export rds file
readr::write_rds(merged_ref_sce, opt$output_reference_file)
