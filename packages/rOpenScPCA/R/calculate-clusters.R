calculate_clusters <- function(
    mat,
    algorithm = "louvain",
    weighting = "jaccard",
    nn = 10,
    resolution = 1, # leiden and louvain
    objective_function = "CPM", # leiden only
    random_seed = 2024,
    cluster_args = list()) {
  set.seed(random_seed)

  # Check input arguments
  stopifnot(
    "The `mat` argument must be a matrix." = class(mat) %in% c("Matrix", "dgCMatrix"),
    "The `mat` matrix must be numeric." = type(mat) %in% c("numeric", "double"),
    "The `mat` matrix must have row names representing cell ids (e.g. barcodes)." = !(is.null(rownames(mat)))
  )

  algorithm <- toLower(algorithm)
  stopifnot(
    "`algorithm` must be one of 'louvain' (default), 'walktrap' or 'leiden'." =
      algorithm %in% c("louvain", "walktrap", "leiden")
  )

  weighting <- toLower(weighting)
  stopifnot(
    "`weighting` must be one of 'jaccard' (default) or 'rank'." =
      weighting %in% c("jaccard", "rank", "leiden")
  )

  # TODO: consider adding more specific checks here?
  stopifnot(
    "`cluster_args` must be a list." = class(cluster_args) == "list"
  )

  # Update cluster_args list with settings that users can directly provide
  if (algorithm != "walktrap") {
    cluster_args$resolution <- resolution
  }
  if (algorithm == "leiden") {
    cluster_args$resolution <- objective_function
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
  # TODO: Should this have _all_ parameters in cluster.args?
  return(
    tibble::tibble(
      barcode = rownames(mat),
      cluster = clusters,
      algorithm = algorithm,
      weighting = weighting,
      nn = nn
    )
  )
}
