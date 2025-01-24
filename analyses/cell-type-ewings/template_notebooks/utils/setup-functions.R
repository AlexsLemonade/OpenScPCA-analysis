# These functions are used in `celltype-exploration.Rmd` 
# They are used for reading in and setting up the cell type results

#' Combine workflow results into a single data frame
#'
#' @param sce Processed SingleCellExperiment object with UMAP embeddings 
#' @param singler_df Data frame with results from `aucell-singler-annotation.sh` workflow
#' @param cluster_df Data frame with results from `evaluate-clusters.sh` workflow
#' @param aucell_df Data frame with results from `run-aucell-ews-signatures.sh` workflow
#' @param cluster_nn Value of nearest neighbors to use for cluster results. Default is 20.
#' @param cluster_res Value of resolution to use for cluster results. Default is 20.  
#'
prep_results <- function(
    sce,
    singler_df,
    cluster_df,
    aucell_df,
    cluster_nn = 20,
    cluster_res = 0.5
) {
  
  ## grab UMAP 
  umap_df <- sce |>
    scuttle::makePerCellDF(use.dimred = "UMAP") |>
    # replace UMAP.1 with UMAP1 and get rid of excess columns
    dplyr::select(barcodes, UMAP1 = UMAP.1, UMAP2 = UMAP.2)
  
  ## prep singler data
  singler_df <- singler_df |> 
    dplyr::mutate(
      # first grab anything that is tumor and label it tumor
      # NA should be unknown
      singler_annotation = dplyr::case_when(
        stringr::str_detect(singler_annotation, "tumor") ~ "tumor",
        is.na(singler_annotation) ~ "unknown", # make sure to separate out unknown labels
        .default = singler_annotation
      ) |>
        forcats::fct_relevel("tumor", after = 0),
      # get the top cell types for plotting later
      singler_lumped = singler_annotation |>
        forcats::fct_lump_n(7, other_level = "All remaining cell types", ties.method = "first") |>
        forcats::fct_infreq() |>
        forcats::fct_relevel("All remaining cell types", after = Inf)
    )
  
  ## prep cluster data 
  cluster_df <- cluster_df |> 
    # filter to the clustering results we want to use 
    dplyr::filter(
      cluster_method == "leiden_mod",
      nn == cluster_nn,
      resolution == cluster_res
    ) |> 
    dplyr::select(
      barcodes = cell_id,
      cluster
    )
  
  ## prep aucell 
  aucell_wide_df <- aucell_df |> 
    dplyr::mutate(
      assignment = auc > auc_threshold
    ) |> 
    tidyr::pivot_wider(
      id_cols = "barcodes",
      names_from = "gene_set",
      values_from = c(auc, assignment)
    )
  
  ## combine into one data frame 
  all_results_df <- umap_df |> 
    dplyr::left_join(singler_df, by = c("barcodes")) |> 
    dplyr::left_join(cluster_df, by = c("barcodes")) |> 
    dplyr::left_join(aucell_wide_df, by = c("barcodes"))
  
  return(all_results_df)
  
}
