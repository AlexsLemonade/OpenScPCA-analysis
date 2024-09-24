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
#' @return Expanded `cluster_df` data frame with these additional columns:
#' - `silhouette_width`, the cell's silhouette width
#' - `other`, the closest cluster other than the one to which the given cell was assigned
#' For more information, see documentation for `bluster::approxSilhouette()`
#'
#' @export
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
    dplyr::rename("silhouette_width" = "width")

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
#' @return Expanded `cluster_df` data frame with these additional columns:
#' - `purity`, the cell's neighborhood purity
#' - `maximum`, the cluster with the highest proportion of observations neighboring the given cell.
#' For more information, see documentation for `bluster::neighborPurity()`
#'
#' @export
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
    tibble::rownames_to_column("cell_id")

  # join with cluster_df in this direction, so that columns in
  # cluster_df come first
  purity_df <- cluster_df |>
    dplyr::inner_join(purity_df, by = c("cell_id"))

  return(purity_df)
}



calculate_stability <- function(
    x,
    clusters,
    algorithm = c("louvain", "walktrap", "leiden"),
    weighting = c("jaccard", "rank", "number"),
    nn = 10,
    resolution = 1,
    objective_function = c("CPM", "modularity"),
    cluster_args = list(),
    n_iter = 20,
    threads = 1,
    seed = NULL,
    pc_name = NULL) {
  # ensure we have a matrix
  pca_matrix <- prepare_pc_matrix(x, pc_name = pc_name)

  # calculate ARI for each cluster result bootstrap replicate
  ari_values <- 1:n_iter |>
    purrr::map(
      \(i) {
        sample_cells <- sample(nrow(pca_matrix), nrow(pca_matrix), replace = TRUE)
        resampled_pca <- pca_matrix[sample_cells, , drop = FALSE]
        original_clusters <- clusters[sample_cells]

        resampled_clusters <- calculate_clusters(
          resampled_pca,
          algorithm,
          weighting,
          nn,
          resolution = resolution,
          objective_function = objective_function,
          cluster_args = cluster_args,
          threads = threads,
          seed = seed
        )

        # return ARI between new clustering and original clustering
        pdfCluster::adj.rand.index(resampled_clusters$cluster, clusters[sample_cells])
      }
    ) |>
    purrr::reduce(c)


  return(
    data.frame(
      replicate = 1:n_iter,
      ari = ari_values
    )
  )
}
