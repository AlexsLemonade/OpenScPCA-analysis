# This file contains functions consumed by scripts in this module.


#' Generate a reference table when testing
#'
#' Creates a data frame with 20% reference cells to run with test data
#' 
#' @param barcodes barcodes in query
#'
#' @returns annotations data frame
generate_test_annotation_df <- function(barcodes){
  n_ref <- floor(length(barcodes) * 0.2)
  n_query <- length(barcodes) - n_ref
  data.frame(
    annotation = c(rep("reference", n_ref), rep("unknown", n_query)),
    row.names = barcodes
  )
}   



#' Prepare and export an internal reference inferCNV annotation file
#'
#' @param reference_name Reference name to create, as recorded in the reference_celltype_tsv
#' @param reference_celltype_tsv TSV with reference names and associated consensus cell types
#' @param celltype_tsv TSV with per-cell cell types to consider when creating the reference
#' @param library_id Library id for the SCE being processed
#' @param testing Logical for whether we are running with test data
prepare_internal_reference_annotations <- function(
    reference_group,
    reference_celltype_tsv,
    celltype_tsv,
    library_id,
    testing) {

  # TODO: currently the code assumes this contains columns `ewing_annotation` and `consensus_annotation`
  celltype_df <- readr::read_tsv(celltype_tsv) |>
    dplyr::mutate(barcodes = glue::glue("{library_id}-{barcodes}"))
  
  # If we're testing, assign 20% of cells to the reference. Otherwise, use cell types categories appropriately
  if (testing) {
    annotation_df <- generate_test_annotation_df(celltype_df$barcodes)
  } else {
    reference_celltypes <- readr::read_tsv(reference_celltype_tsv) |>
      dplyr::filter(reference_name == reference_group) |>
      dplyr::pull(consensus_celltype)
    stopifnot("Could not find any cell types to include in the specified internal reference." = 
                length(reference_celltypes) > 0)

    annotation_df <- data.frame(
      annotation = ifelse(
        # TODO: If we want to use non-ewings projects in this script, this will NOT WORK.
        # We may need project-specific helper functions for this spot in the code if/when the time comes.
        consensus_annotation %in% reference_celltypes & !stringr::str_detect(ewing_annotation, "tumor"),
        "reference",
        "unknown"
      ),
      row.names = celltype_df$barcodes
    )
  }
  annotation_df
}


# # Export annotations file
# prepare_pooled_reference_annotations(
#   colnames(raw_counts_matrix), # all cell ids in the matrix
#   colnames(ref_sce), # reference names specifically
#   opts$annotation_file, # output file
#   opts$testing # logical if we're running with test data
# )

#' Prepare and export a pooled reference inferCNV annotation file
#'
#' @param all_cell_ids All cell ids present in the inferCNV input
#' @param reference_cell_ids Cell ids present in the reference set of cells
#' @param testing Logical for whether we are running with test data
prepare_pooled_reference_annotations <- function(
    all_cell_ids,
    reference_cell_ids,
    annotation_file, 
    testing) {
  if (testing) {
    annotation_df <- generate_test_annotation_df(all_cell_ids)
  } else {
    # "unknown" cells are uncharacterized, and "reference" cells are in the reference
    annotation_df <- data.frame(
      annotations = dplyr::if_else(
        all_cell_ids %in% reference_cell_ids, "reference", "unknown"
      ), 
      row.names = all_cell_ids
    )
  }
  annotation_df
}


# Helper function to remove unneeded slots
# from a reference SCE to save space
clean_sce <- function(sce) {
  logcounts(sce) <- NULL
  assay(sce, "spliced") <- NULL
  reducedDim(sce, "PCA") <- NULL
  reducedDim(sce, "UMAP") <- NULL
  
  # ensure the counts matrix is sparse
  counts(sce) <- as(counts(sce), "CsparseMatrix")
  
  return(sce)
}

# Helper function to join consensus cell types into a colData slot
# - sce_coldata: the colData slot of an SCE object, which will be updated and returned
# - consensus_celltype_df: a data frame with at least columns `sce_cell_id` and `consensus_annotation`
consensus_to_coldata <- function(sce_coldata, consensus_celltype_df) {
  sce_coldata |>
    as.data.frame() |>
    # temporarily make the rownames a column so we can join consensus
    tibble::rownames_to_column(var = "sce_cell_id") |>
    dplyr::left_join(consensus_celltype_df, by = "sce_cell_id") |>
    dplyr::select(-sce_cell_id) |>
    # make it a DataFrame again
    DataFrame(row.names = rownames(sce_coldata))
}
