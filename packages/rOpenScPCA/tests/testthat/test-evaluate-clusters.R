suppressPackageStartupMessages(library(SingleCellExperiment))

set.seed(2024)
sce <- splatter::simpleSimulate(nGenes = 1000, verbose = FALSE) |>
  scater::logNormCounts() |>
  scater::runPCA(ncomponents = 10)
test_mat <- reducedDim(sce, "PCA")


cluster_df <- calculate_clusters(test_mat)

test_that("calculate_silhouette works as expected", {
  df <- calculate_silhouette(test_mat, cluster_df)

  expect_setequal(
    colnames(df),
    c(colnames(cluster_df), "silhouette_width", "other")
  )
  expect_equal(df$cell_id, rownames(test_mat))
  expect_equal(df$cluster, cluster_df$cluster)
  expect_vector(df$silhouette_width, ptype = numeric())
  expect_s3_class(df$other, "factor")
})



test_that("calculate_purity works as expected", {
  df <- calculate_purity(test_mat, cluster_df)

  expect_setequal(
    colnames(df),
    c(colnames(cluster_df), "purity", "maximum")
  )
  expect_equal(df$cell_id, rownames(test_mat))
  expect_equal(df$cluster, cluster_df$cluster)
  expect_vector(df$purity, ptype = numeric())
  expect_s3_class(df$maximum, "factor")
})





test_that("calculate_stability works as expected with defaults", {
  # note that we suppress warnings since this calculation done on fake
  # test data gives expected warnings about ties during the ARI calculation.
  suppressWarnings({
    df <- calculate_stability(test_mat, cluster_df$cluster)
  })

  expected_names <- colnames(cluster_df)[!(colnames(cluster_df) %in% c("cell_id", "cluster"))]
  expect_setequal(
    colnames(df),
    c("replicate", "ari", expected_names)
  )
  expect_equal(df$replicate, 1:20) # checks rows too
  expect_vector(df$ari, ptype = numeric())
})


test_that("calculate_stability works as expected with different replicates", {
  # note that we suppress warnings since this calculation done on fake
  # test data gives expected warnings about ties during the ARI calculation.
  suppressWarnings({
    df <- calculate_stability(test_mat, cluster_df$cluster, replicates = 2)
  })
  expect_equal(nrow(df), 2)
})



test_that("calculate_stability works as expected with object and pc_name", {
  reducedDimNames(sce) <- "my_pca"

  # note that we suppress warnings since this calculation done on fake
  # test data gives expected warnings about ties during the ARI calculation.
  suppressWarnings({
    df <- calculate_stability(
      sce,
      cluster_df$cluster,
      replicates = 2,
      pc_name = "my_pca"
    )
  })
  expect_equal(nrow(df), 2)
})



test_that("calculate_stability errors as expected", {
  expect_error({
    # mismatched cluster vector length
    calculate_stability(test_mat, cluster_df$cluster[1:5])
  })

  expect_error({
    # cluster_df not a vector
    calculate_stability(test_mat, cluster_df)
  })


  # the next test ensures we fail even if the length of clusters
  # is the same as the nrow in test_mat

  # make it smaller for test
  test_mat_small <- test_mat[1:5, ]

  # make df with same length as rows in test_mat
  bad_df <- data.frame(
    cluster = sample(1:2, nrow(test_mat_small), replace = T),
    cell_id = rownames(test_mat_small),
    col1 = 1,
    col2 = 1,
    col3 = 1
  )

  expect_error({
    # cluster_df not a vector, even though df length matches test mat
    calculate_stability(test_mat_small, bad_df)
  })
})
