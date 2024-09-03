test_that("calculate_clusters runs with defaults", {
  test_mat <- matrix(
    runif(1000, 0.1, 10),
    nrow = 100,
    ncol = 10
  )

  rownames(test_mat) <- 1:100
  cluster_df <- calculate_clusters(test_mat)


  expect_equal(
    names(cluster_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn")
  )
  expect_equal(
    cluster_df$cell_id,
    as.character(1:100)
  )
  expect_true(
    is.factor(cluster_df$cluster)
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
})
