#' Calculate graph-based clusters from a provided matrix
#'
#' @param mat Matrix, usually of PCs, where each row is a cell. Matrix must have rownames of cell ids (e.g., barcodes)
#' @param algorithm Clustering algorithm to use. Must be one of "louvain" (default), "walktrap", or "leiden"
#' @param weighting Weighting scheme to use. Must be one of "jaccard" (default) or "rank"
#' @param nn Number of nearest neighbors. Default is 10.
#' @param resolution Resolution parameter used by louvain and leiden clustering only. Default is 1.
#' @param objective_function Leiden-specific parameter for whether to use the Constant Potts Model ("CPM"; default) or "modularity"
#' @param seed Random seed to set for clustering. Default is 2024.
#' @param cluster_args List of additional arguments to pass to the chosen clustering function. See igraph documentation for details: https://igraph.org/r/html/latest
#'
#' @return A data frame of cluster results with columns `cell_id` and `cluster`. Additional columns represent algorithm parameters
#'   and include at least: `algorithm`, `weighting`, and `nn`. Louvain and leiden clustering will also include `resolution`, and
#'   leiden clustering will further include `objective_function`. Any other values passed to `cluster_args` will also be included.
#' @export
#'
#' @examples
#' \dontrun{
#' # cluster using default parameters
#' cluster_df <- calculate_clusters(pca_matrix)
#'
#' # cluster using the leiden algorithm with a resolution of 0.1
#' cluster_df <- calculate_clusters(pca_matrix, algorithm = "leiden", resolution = 0.1)
#'
#' # cluster using the leiden algorithm with a non-default of 3 iterations
#' cluster_df <- calculate_clusters(
#'   pca_matrix,
#'   algorithm = "leiden",
#'   cluster_args = list(n_iterations = 3)
#' )
#' }
calculate_clusters <- function(
    mat,
    algorithm = c("louvain", "walktrap", "leiden"),
    weighting = c("jaccard", "rank"),
    nn = 10,
    resolution = 1, # louvain or leiden
    objective_function = c("CPM", "modularity"), # leiden only
    cluster_args = list(),
    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check input arguments
  stopifnot(
    "The `mat` argument must be a matrix." = any(class(mat) %in% c("matrix", "Matrix")),
    "The `mat` matrix must have row names representing cell ids, e.g. barcodes." = !(is.null(rownames(mat)))
  )

  algorithm <- match.arg(algorithm)
  weighting <- match.arg(weighting)
  objective_function <- match.arg(objective_function)

  stopifnot(
    "`resolution` must be numeric" = is.numeric(resolution),
    "`nn` must be numeric" = is.numeric(nn)
  )

  if (length(cluster_args) != 0) {
    stopifnot(
      "`cluster_args` must be a named list." = is.list(cluster_args) && !(is.null(names(cluster_args)))
    )
  }

  # Update cluster_args list with parameters that users can directly provide
  # note that clusterRows throws an error if this list has a param not used by the chosen algorithm
  if (algorithm != "walktrap") {
    cluster_args$resolution <- resolution
  }
  if (algorithm == "leiden") {
    cluster_args$objective_function <- objective_function
  }


  # Perform clustering
  clusters <- bluster::clusterRows(
    mat,
    bluster::NNGraphParam(
      k = nn,
      type = weighting,
      cluster.fun = algorithm,
      cluster.args = cluster_args
    )
  )


  # Transform results into a table and return
  cluster_df <- data.frame(
    cell_id = rownames(mat),
    cluster = clusters,
    algorithm = algorithm,
    weighting = weighting,
    nn = nn
  ) |>
    dplyr::bind_cols(
      dplyr::bind_rows(cluster_args)
    )

  return(cluster_df)
}
