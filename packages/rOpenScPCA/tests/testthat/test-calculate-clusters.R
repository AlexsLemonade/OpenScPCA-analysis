test_mat <- matrix(
  runif(1000, 0.1, 10),
  nrow = 100,
  ncol = 10
)
rownames(test_mat) <- as.character(1:100)

test_that("calculate_clusters runs with defaults", {
  cluster_df <- calculate_clusters(test_mat)

  expect_equal(
    names(cluster_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
  )
  expect_equal(
    cluster_df$cell_id,
    as.character(1:100)
  )
  expect_true(
    is.factor(cluster_df$cluster)
  )
  expect_true(
    all(cluster_df$algorithm == "louvain")
  )
  expect_true(
    all(cluster_df$weighting == "jaccard")
  )
  expect_true(
    all(cluster_df$nn == 10)
  )
  expect_true(
    all(cluster_df$resolution == 1)
  )
})


test_that("calculate_clusters runs with additional cluster_args", {
  cluster_df <- calculate_clusters(
    test_mat,
    algorithm = "leiden",
    # the default is 2
    cluster_args = list(n_iterations = 3)
  )

  expect_setequal(
    names(cluster_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution", "objective_function", "n_iterations")
  )
  expect_true(
    all(cluster_df$n_iterations == 3)
  )
})




test_that("calculate_clusters errors as expected", {
  test_mat_nonames <- test_mat
  rownames(test_mat_nonames) <- NULL

  expect_error(
    calculate_clusters(test_mat_nonames),
    "The `mat` matrix must have row names representing cell ids, e.g. barcodes."
  )
  expect_error(
    calculate_clusters("not a matrix"),
    "The `mat` argument must be a matrix."
  )

  expect_error(
    calculate_clusters(test_mat, resolution = "string"),
    "`resolution` must be numeric"
  )
  expect_error(
    calculate_clusters(test_mat, nn = "string"),
    "`nn` must be numeric"
  )
})
