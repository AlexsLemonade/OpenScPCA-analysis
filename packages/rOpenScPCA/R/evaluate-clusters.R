#' Calculate the silhouette width of clusters
#'
#' This function uses the `bluster::approxSilhouette()` function to calculate the
#' silhouette width for a clustering result. These results can be used downstream to
#' calculate the average silhouette width, a popular metric for cluster evaluation.
#'
#' @param x Either a matrix of principal components (PCs), or a SingleCellExperiment
#'  or Seurat object containing PCs. If a matrix is provided, rows should be cells
#'  and columns should be PCs, and row names should be cell ids (e.g., barcodes).
#' @param cluster_df A data frame that contains at least the columns `cell_id` and
#'  `cluster`. The `cell_id` values should match either the PC matrix row names,
#'  or the SingleCellExperiment/Seurat object cell ids. Typically this will be output from
#'  the `rOpenScPCA::calculate_clusters()` function.
#' @param pc_name Optionally, the name of the PC matrix in the object. Not used if a
#'   matrix is provided. If the name is not provided, the name "PCA" is assumed for
#'   SingleCellExperiment objects, and "pca" for Seurat objects.
#'
#' @return Expanded `cluster_df` data frame with the additional column `silhouette_width`
#' @export
#'
#' @examples
#' \dontrun{
#' # calculate silhouette width for clusters stored in a data frame
#' cluster_df <- calculate_silhouette(sce_object, cluster_df)
#' }
calculate_silhouette <- function(
    x,
    cluster_df,
    pc_name = NULL) {
  x <- prepare_pc_matrix(x, pc_name)

  expected_df_names <- c("cell_id", "cluster")
  stopifnot(
    "Expected columns 'cell_id' and 'cluster' in the cluster_df." =
      all(expected_df_names %in% colnames(cluster_df))
  )

  silhouette_df <- x |>
    bluster::approxSilhouette(cluster_df$cluster) |>
    as.data.frame() |>
    tibble::rownames_to_column("cell_id") |>
    dplyr::select(
      cell_id,
      cluster,
      # TODO: should we keep `other`? If so, this can just be
      # dplyr::rename instead of select
      silhouette_width = width
    )

  # join with cluster_df in this direction, so that columns in
  # cluster_df come first
  silhouette_df <- cluster_df |>
    dplyr::inner_join(silhouette_df, by = c("cell_id", "cluster"))

  return(silhouette_df)
}




#' Calculate the neighborhood purity of clusters
#'
#' This function uses the `bluster::neighborPurity()` function to calculate the
#' neighborhood purity values for a clustering result.
#'
#' @param x Either a matrix of principal components (PCs), or a SingleCellExperiment
#'  or Seurat object containing PCs. If a matrix is provided, rows should be cells
#'  and columns should be PCs, and row names should be cell ids (e.g., barcodes).
#' @param cluster_df A data frame that contains at least the columns `cell_id` and
#'  `cluster`. The `cell_id` values should match either the PC matrix row names,
#'  or the SingleCellExperiment/Seurat object cell ids. Typically this will be output from
#'  the `rOpenScPCA::calculate_clusters()` function.
#' @param pc_name Optionally, the name of the PC matrix in the object. Not used if a
#'   matrix is provided. If the name is not provided, the name "PCA" is assumed for
#'   SingleCellExperiment objects, and "pca" for Seurat objects.
#' @param ... Additional arguments to pass to `bluster::neighborPurity()`
#'
#' @return Expanded `cluster_df` data frame with the additional column `purity`
#' @export
#'
#' @examples
#' \dontrun{
#' # calculate neighborhood purity for clusters stored in a data frame
#' cluster_df <- calculate_purity(sce_object, cluster_df)
#' }
calculate_purity <- function(
    x,
    cluster_df,
    pc_name = NULL,
    ...) {
  x <- prepare_pc_matrix(x, pc_name)

  expected_df_names <- c("cell_id", "cluster")
  stopifnot(
    "Expected columns 'cell_id' and 'cluster' in cluster_df." =
      all(expected_df_names %in% colnames(cluster_df))
  )

  purity_df <- x |>
    bluster::neighborPurity(cluster_df$cluster) |>
    as.data.frame() |>
    tibble::rownames_to_column("cell_id") |>
    # TODO: should we keep `maximum`? If so, can remove this select statement
    dplyr::select(-maximum)

  # join with cluster_df in this direction, so that columns in
  # cluster_df come first
  purity_df <- cluster_df |>
    dplyr::inner_join(purity_df, by = c("cell_id"))

  return(purity_df)
}
