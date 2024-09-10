test_mat <- matrix(
  runif(1000, -3, 3),
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

  expect_s3_class(
    cluster_df$cluster,
    "factor"
  )

  expect_equal(
    unique(cluster_df$algorithm),
    "louvain"
  )
  expect_equal(
    unique(cluster_df$weighting),
    "jaccard"
  )
  expect_equal(
    unique(cluster_df$nn),
    10
  )
  expect_equal(
    unique(cluster_df$resolution),
    1
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
  expect_equal(
    unique(cluster_df$n_iterations),
    3
  )
})




test_that("calculate_clusters errors as expected", {
  test_mat_nonames <- test_mat
  rownames(test_mat_nonames) <- NULL

  expect_error(calculate_clusters(test_mat_nonames))
  expect_error(calculate_clusters("not a matrix"))
  expect_error(calculate_clusters(test_mat, resolution = "string"))
  expect_error(calculate_clusters(test_mat, nn = "string"))
})
