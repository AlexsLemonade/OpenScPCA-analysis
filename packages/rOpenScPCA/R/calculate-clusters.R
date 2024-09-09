#' Calculate graph-based clusters from a provided matrix
#'
#' @param mat Matrix, usually of PCs, where each row is a cell. Matrix must have rownames of cell ids (e.g., barcodes)
#' @param algorithm Clustering algorithm to use. Must be one of "louvain" (default), "walktrap", or "leiden"
#' @param weighting Weighting scheme to use. Must be one of "jaccard" (default) or "rank"
#' @param nn Number of nearest neighbors. Default is 10.
#' @param resolution Resolution parameter used by louvain and leiden clustering only. Default is 1.
#' @param objective_function Leiden-specific parameter for whether to use the Constant Potts Model ("CPM"; default) or "modularity"
#' @param seed Random seed to set for clustering. Default is 2024.
#' @param cluster_args List of additional arguments to pass to the clustering function.
#'
#' @return A data frame of cluster results with columns `cell_id` and `cluster`. Additional columns (TODO!) represent algorithm parameters.
#' @export
#'
#' @examples
#' \dontrun{
#' # TODO WILL ADD EXAMPLES ONCE ARGS ARE FINALIZED
#' }
calculate_clusters <- function(
    mat,
    algorithm = c("louvain", "walktrap", "leiden"),
    weighting = c("jaccard", "rank"),
    nn = 10,
    resolution = 1, # louvain or leiden
    objective_function = "CPM", # leiden only
    cluster_args = list(),
    seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Check input arguments
  stopifnot(
    "The `mat` argument must be a matrix." = any(class(mat) %in% c("matrix", "Matrix")),
    "The `mat` matrix must have row names representing cell ids (e.g. barcodes)." = !(is.null(rownames(mat)))
  )

  algorithm <- match.arg(algorithm)
  print(algorithm)
  stopifnot(
    "`algorithm` must be one of 'louvain' (default), 'walktrap' or 'leiden'." =
      algorithm %in% c("louvain", "walktrap", "leiden")
  )

  weighting <- match.arg(weighting)
  stopifnot(
    "`weighting` must be one of 'jaccard' (default) or 'rank'." =
      weighting %in% c("jaccard", "rank")
  )

  # TODO: consider adding more specific checks here?
  stopifnot(
    "`cluster_args` must be a list." = is.list(cluster_args)
  )

  # Update cluster_args list with settings that users can directly provide
  # clusterRows throws an error if this list has a param not used by the chosen algorithm
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
  # TODO: Should this have _all_ (non-default/user-specified) parameters in cluster.args?
  return(
    data.frame(
      cell_id = rownames(mat),
      cluster = clusters,
      algorithm = algorithm,
      weighting = weighting,
      nn = nn
    )
  )
}
