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



#' Calculate cluster stability using the Adjusted Rand Index (ARI)
#'
#' This function generates and clusters, using provided parameters, bootstrap
#' replicates calculates the Adjusted Rand Index (ARI) between each set of bootstrapped
#' clusters and the original provided clusters. ARI measures similarity between different
#' cluster results, where a value of 0 indicates an entirely random relationship between
#' results, and a value of 1 indicates perfect agreement.
#'
#' When assessing stability, you should specify the same clustering parameters here as
#' were used to calculate the original clusters.
#'
#' Note that this function will also make use of bluster::clusterRows() with the
#' bluster::NNGraphParam() function on a principal components matrix. Note that defaults
#' for some arguments may differ from the bluster::NNGraphParam() defaults.
#' Specifically, the clustering algorithm defaults to "louvain" and the weighting scheme
#' to "jaccard" to align with common practice in scRNA-seq analysis.
#'
#'
#' @param x An object containing PCs that clusters were calculated from. This can be
#'   either a SingleCellExperiment object, a Seurat object, or a matrix where columns
#'   are PCs and rows are cells. If a matrix is provided, it must have row names of cell
#'   ids (e.g., barcodes).
#' @param clusters A vector of cluster ids, typically a numeric factor variable, obtained
#'   by previously clustering the PCs.
#' @param replicates Number of bootstrap replicates to perform. Default is 20.
#' @param seed Random seed
#' @param pc_name Optionally, the name of the PC matrix in the object. Not used if a
#'   matrix is provided. If the name is not provided, the name "PCA" is assumed for
#'   SingleCellExperiment objects, and "pca" for Seurat objects.
#' @param ... Additional arguments to pass to `calculate_clusters()` which calculates
#'   bootstrapped clusters. Usually, these will be the same arguments used to generate
#'   the original clusters.
#'
#' @return Data frame with columns `replicate` and `ari`, representing the given bootstrap replicate
#'   and its ARI value, respectively, and columns representing clustering algorithm parameters which
#'   include at least `algorithm`, `weighting`, and `nn`. Louvain and leiden clustering will also
#'   include `resolution`, and leiden clustering will further include `objective_function`.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # First, cluster PCs from a SingleCellExperiment object using default parameters
#' # and setting a seed for reproducibility
#' cluster_df <- calculate_clusters(sce_object, seed = 11)
#' # Second, calculate cluster stability using default parameters
#' stability_df <- calculate_stability(sce_object, cluster_df$clusters, seed = 11)
#'
#'
#' # First, cluster PCs from a SingleCellExperiment object using default parameters
#' # and setting a seed for reproducibility
#' cluster_df <- calculate_clusters(sce_object, seed = 11)
#' # Second, calculate cluster stability using default parameters and 50 replicates
#' stability_df <- calculate_stability(
#'   sce_object,
#'   cluster_df$clusters,
#'   replicates = 50,
#'   seed = 11
#' )
#'
#'
#' # First, cluster PCs from a SingleCellExperiment object using the leiden
#' # algorithm and a resolution of 0.1
#' cluster_df <- calculate_clusters(
#'   sce_object,
#'   algorithm = "leiden",
#'   resolution = 0.1,
#'   seed = 11
#' )
#' # Second, calculate cluster stability using the same parameters as were used
#' # for the initial clustering
#' stability_df <- calculate_stability(
#'   sce_object,
#'   cluster_df$clusters,
#'   algorithm = "leiden",
#'   resolution = 0.1,
#'   seed = 11
#' )
#' }
calculate_stability <- function(
    x,
    clusters,
    replicates = 20,
    seed = NULL,
    pc_name = NULL,
    ...) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # ensure we have a matrix
  pca_matrix <- prepare_pc_matrix(x, pc_name = pc_name)

  # check clusters and matrix compatibility
  stopifnot(
    "The number of rows in the matrix must equal the length of the clusters vector." =
      nrow(pca_matrix) == length(clusters)
  )

  # calculate ARI for each cluster result bootstrap replicate
  all_ari_df <- 1:replicates |>
    purrr::map(
      \(i) {
        sample_cells <- sample(nrow(pca_matrix), replace = TRUE)
        resampled_pca <- pca_matrix[sample_cells, ]
        original_clusters <- clusters[sample_cells]

        resampled_df <- calculate_clusters(resampled_pca, ...)

        ari <- pdfCluster::adj.rand.index(resampled_df$cluster, original_clusters)

        # return df with ari and clustering parameters
        ari_df <- resampled_df |>
          dplyr::slice(1) |>
          dplyr::select(!c("cell_id", "cluster")) |>
          dplyr::mutate(
            # define this variable here to ensure it's numeric
            replicate = i,
            ari = ari,
            # ensure these columns come first
            .before = "algorithm"
          )

        return(ari_df)
      }
    ) |>
    dplyr::bind_rows()


  return(all_ari_df)
}
