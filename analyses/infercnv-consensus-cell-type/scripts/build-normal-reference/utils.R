# Utility functions for building normal references

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
