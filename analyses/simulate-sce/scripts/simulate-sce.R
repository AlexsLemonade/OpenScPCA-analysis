#!/usr/bin/env Rscript

# Script to generate simulated single-cell data from a real dataset

# Works on a single sample directory at a time, processing all SCE files within that directory.


# Parse arguments --------------------------------------------------------------

library(optparse)

# Define and parse the command line arguments
option_list <- list(
  make_option(
    c("-s", "--sample_dir"),
    type = "character",
    default = NULL,
    help = "The sample directory for the real data to use for simulation"
  ),
  make_option(
    c("-m", "--metadata_file"),
    type = "character",
    default = NULL,
    help = "The permuted metadata file to use for metadata replacement"
  ),
  make_option(
    c("-o", "--output_dir"),
    type = "character",
    default = NULL,
    help = "The output directory. Output files will be given the same names as input"
  ),
  make_option(
    c("-n", "--ncells"),
    type = "integer",
    default = 100,
    help = "The number of cells to simulate."
  ),
  make_option(
    c("--seed"),
    type = "integer",
    default = 2024,
    help = "Random number seed for simulation."
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# Load libraries ---------------------------------------------------------------

# wait to load these until after the command line arguments are parsed
# Some take a while to load, and we'd rather have the script be responsive

# There are some annoying warning here, let's suppress them
suppressPackageStartupMessages({
  suppressWarnings(library(SingleCellExperiment))
})


# Functions --------------------------------------------------------------------

#' Randomize labels
#'
#' A function to randomly select labels from a set,
#' while ensuring each label is included at least once, if possible.
#'
#' @param label_set a vector a labels
#' @param n the number of labels to select
#'
#' @return a vector of length n
#'
#' @examples
random_label <- function(label_set, n) {
  # ensure labels are unique
  label_set <- unique(label_set)
  if (length(label_set) > n) {
    warning("The number of labels is greater than the number requested; not all labels will be used.")
  }
  r_labels <- sample(label_set, n, replace = TRUE)
  # add each label at least once
  idx <- sample.int(n, min(n, length(label_set)))
  r_labels[idx] <- label_set
  return(r_labels)
}

#' Simulate an SCE object based on a previous object, retaining metadata
#'
#' @param sce a SingleCellExperiment object, with the formatted expected from ScPCA
#' @param ncells The number of cells to simulate
#' @param replacement_metadata A data frame with sample metadata to replace in the sce object
#' @param processed Boolean indicating whether the data is processed or not

#'
#' @return a simulated SingleCellExperiment object with metadata matching the input
#' @export
#'
#' @import SingleCellExperiment
#'
#' @examples
simulate_sce <- function(sce, ncells, replacement_metadata, processed) {
  # check parameters
  stopifnot(
    "sce must be a SingleCellExperiment" = is(sce, "SingleCellExperiment"),
    "ncells must be a positive integer" = is.integer(ncells) && ncells > 0,
    "replacement_metadata should be a data frame with the same number of rows as the sce sample_metadata" = (
      is.data.frame(replacement_metadata) && nrow(replacement_metadata) == nrow(metadata(sce)$sample_metadata)
    ),
    "processed must be a boolean" = is.logical(processed)
  )

  # define a subset of cells to simulate
  ncells <- min(ncells, ncol(sce))
  cell_subset <- sample(colnames(sce), ncells)
  sce_sim <- sce[, cell_subset]

  ### Reduce and remove metadata -----------------------------------------------

  # match column names from permuted metadata to sample_metadata
  replacement_metadata <- replacement_metadata |>
    dplyr::rename(
      sample_id = scpca_sample_id
    ) |>
    dplyr::select(any_of(colnames(metadata(sce)$sample_metadata))) |>
    dplyr::mutate( # convert age to the type in the table before replacement
      age = as(age, type(metadata(sce)$sample_metadata$age))
    )

  # replace sample metadata fields with permuted values
  metadata(sce_sim)$sample_metadata <- metadata(sce_sim)$sample_metadata |>
    dplyr::rows_update(replacement_metadata, by = "sample_id")

  # remove miQC model that may exist
  metadata(sce_sim)$miQC_model <- NULL

  # reduce the cell type data matrices, if present
  if (!is.null(metadata(sce_sim)$singler_results)) {
    sim_cells <- cell_subset[cell_subset %in% rownames(metadata(sce_sim)$singler_results)]
    singler_results <- metadata(sce_sim)$singler_results[sim_cells, ]
    # remove any gene data present to save space
    metadata(singler_results)$de.genes <- NULL
    metadata(singler_results)$common.genes <- NULL
    metadata(sce_sim)$singler_results <- singler_results

  }
  if (!is.null(metadata(sce_sim)$cellassign_predictions)) {
    cellassign_subset <- metadata(sce_sim)$cellassign_predictions |>
      dplyr::filter(barcode %in% cell_subset)
    metadata(sce_sim)$cellassign_predictions <- cellassign_subset
  }

  # Adjust cluster/cell type labels --------------------------------------------
  if (!is.null(colData(sce)$submitter_celltype_annotation)) {
    submitter_set <- unique(colData(sce)$submitter_celltype_annotation)
    colData(sce_sim)$submitter_celltype_annotation <- random_label(submitter_set, ncells)
  }
  # clusters
  if (!is.null(colData(sce)$cluster)) {
    cluster_set <- unique(colData(sce)$cluster)
    colData(sce_sim)$cluster <- random_label(cluster_set, ncells)
  }
  # cellassign cell types
  if (!is.null(colData(sce)$cellassign_celltype_annotation)) {
    cellassign_set <- unique(colData(sce)$cellassign_celltype_annotation)
    colData(sce_sim)$cellassign_celltype_annotation <- random_label(cellassign_set, ncells)
  }
  # singler cell types
  if (!is.null(colData(sce)$singler_celltype_ontology)) {
    # create a mapping for singler ontology and annotation
    singler_df <- colData(sce) |>
      data.frame() |>
      dplyr::select(singler_celltype_ontology, singler_celltype_annotation) |>
      dplyr::distinct()
    singler_dict <- setNames(
      singler_df$singler_celltype_annotation,
      singler_df$singler_celltype_ontology
    )
    colData(sce_sim)$singler_celltype_ontology <- random_label(names(singler_dict), ncells)
    # use dictionary to get proper annotations
    colData(sce_sim)$singler_celltype_annotation <- unname(
      singler_dict[colData(sce_sim)$singler_celltype_ontology]
    )
  }


  # Perform simulation ---------------------------------------------------------
  # remove any all-zero droplets that might have slipped through
  droplets <- colnames(sce)[which(colSums(counts(sce)) > 0)]
  # use a large subset for estimating parameters, but not all
  est_matrix <- counts(sce)[, sample(droplets, min(1000, ncol(sce)))]
  sim_params <- splatter::simpleEstimate(as.matrix(est_matrix))
  sim_params@nCells <- ncells
  # get spliced ratio
  spliced_ratio <- sum(assay(sce, "spliced")) / sum(counts(sce))
  counts(sce_sim, withDimnames = FALSE) <- counts(
    splatter::simpleSimulate(sim_params, verbose = FALSE)
  )
  # make a spliced assay
  assay(sce_sim, "spliced") <- round(counts(sce_sim) * spliced_ratio)

  # recalculate dimension reduction for processed data
  if (processed) {
    logcounts(sce_sim, withDimnames = FALSE) <- log1p(counts(sce_sim)) # use log1p for speed
    sce_sim <- scater::runPCA(sce_sim, 10) # we don't need all the PCA components
    if ("UMAP" %in% reducedDimNames(sce)) { # only run UMAP if it was run before
      sce_sim <- scater::runUMAP(sce_sim)
    }
  }

  # Add any altExps ------------------------------------------------------------
  altExpNames(sce_sim) |> purrr::walk(\(name) {
    alt_sce <- altExp(sce_sim, name)
    sim_counts <- apply(counts(alt_sce), 1, sample) |> t() # Randomize each row
    counts(alt_sce, withDimnames = FALSE) <- sim_counts
    if (processed) {
      logcounts(alt_sce, withDimnames = FALSE) <- log1p(counts(alt_sce))
    }
    # recalculate rowstats
    rowData(alt_sce)[, c("sum", "detected")] <- scuttle::perFeatureQCMetrics(alt_sce)
    altExp(sce_sim, name) <- alt_sce
  })

  # Replace and update column stats --------------------------------------------
  # store column names to restore order later
  coldata_names <- names(colData(sce_sim))
  rowdata_names <- names(rowData(sce_sim))


  # remove stats that will be recalculated
  remove_colstats <- c(
    "sum",
    "detected",
    "total",
    "subsets_mito_sum",
    names(colData(sce_sim))[grep("altexps_.*_(sum|detected|percent)", names(colData(sce_sim)))]
  )
  colData(sce_sim)[, remove_colstats] <- NULL

  remove_rowstats <- c(
    "sum",
    "detected"
  )
  rowData(sce_sim)[, remove_rowstats] <- NULL

  sce_sim <- sce_sim |>
    scuttle::addPerCellQCMetrics() |>
    scuttle::addPerFeatureQCMetrics()
  colData(sce_sim)$subsets_mito_sum <- colData(sce_sim)$sum * colData(sce_sim)$subsets_mito_percent / 100

  # restore column order
  colData(sce_sim) <- colData(sce_sim)[, coldata_names]
  rowData(sce_sim) <- rowData(sce_sim)[, rowdata_names]

  return(sce_sim)
}


# Main Script body -------------------------------------------------------------

set.seed(opts$seed)

# make sure sex is read as a character to prevent all females -> logical false
metadata <- readr::read_tsv(opts$metadata_file, col_types = readr::cols(sex = "c"))

# get file list
sce_files <- list.files(
  opts$sample_dir,
  pattern = "_(processed|filtered|unfiltered)\\.rds$",
  full.names = TRUE
)

fs::dir_create(opts$output_dir)

# perform simulations for each file
purrr::walk(sce_files, \(sce_file) {
  is_processed <- grepl("_processed\\.rds$", sce_file)
  # load the real data
  real_sce <- readr::read_rds(sce_file)
  # get the matching library metadata for replacing
  replacement_metadata <- metadata |>
    dplyr::filter(scpca_library_id == metadata(real_sce)$library_id)
  # simulate the data
  sim_sce <- simulate_sce(real_sce, opts$ncells, replacement_metadata, processed = is_processed)
  # save the simulated data
  sim_file <- file.path(
    opts$output_dir,
    basename(sce_file)
  )
  readr::write_rds(sim_sce, sim_file, compress = "bz2")
})
