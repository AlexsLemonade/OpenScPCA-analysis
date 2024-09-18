suppressPackageStartupMessages(library(SingleCellExperiment))

set.seed(2024)
sce <- splatter::simpleSimulate(nGenes = 1000, verbose = FALSE) |>
  scater::logNormCounts() |>
  scater::runPCA(ncomponents = 10)

test_mat <- reducedDim(sce, "PCA")


test_that("calculate_clusters_sweep works as expected", {
  cluster_sweep_df <- calculate_clusters_sweep(
    sce,
    nn = c(10, 15),
    resolution = c(0.5, 1)
  )

  expect_setequal(
    colnames(cluster_sweep_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution", "cluster_set")
  )

  expect_equal(
    cluster_sweep_df$cell_id,
    rep(rownames(test_mat), 4)
  )

  expect_s3_class(
    cluster_sweep_df$cluster,
    "factor"
  )

  expect_equal(
    unique(cluster_sweep_df$algorithm),
    "louvain"
  )
  expect_equal(
    unique(cluster_sweep_df$weighting),
    "jaccard"
  )
  expect_setequal(
    unique(cluster_sweep_df$nn),
    c(10, 15)
  )
  expect_setequal(
    unique(cluster_sweep_df$resolution),
    c(0.5, 1)
  )
})



test_that("calculate_clusters_sweep works as expected with non-default algorithm", {
  cluster_sweep_df <- calculate_clusters_sweep(
    sce,
    algorithm = "leiden",
    objective_function = "modularity",
    resolution = c(0.5, 1)
  )

  expect_setequal(
    colnames(cluster_sweep_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution", "objective_function", "cluster_set")
  )


  expect_equal(
    unique(cluster_sweep_df$algorithm),
    "leiden"
  )

  expect_equal(
    unique(cluster_sweep_df$objective_function),
    "modularity"
  )

  expect_setequal(
    unique(cluster_sweep_df$resolution),
    c(0.5, 1)
  )
})
