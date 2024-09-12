suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(Seurat))

sce <- scpcaTools:::sim_sce(n_cells = 1000, n_genes = 50, n_empty = 0)
test_mat <- matrix(
  runif(1000, -3, 3),
  nrow = 1000,
  ncol = 50
)
rownames(test_mat) <- colnames(sce)
colnames(test_mat) <- paste0("PC_", 1:50) # quiet seurat warnings
reducedDim(sce, "PCA") <- test_mat

srat <- CreateSeuratObject(counts = counts(sce), assay = "RNA")
srat[["pca"]] <- CreateDimReducObject(
  embeddings = test_mat,
  key = "PC_",
  assay = "RNA"
)



test_that("calculate_clusters runs with defaults", {
  cluster_df <- calculate_clusters(test_mat)

  expect_equal(
    names(cluster_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
  )
  expect_equal(
    cluster_df$cell_id,
    rownames(test_mat)
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
  expect_error(
    calculate_clusters(
      test_mat,
      cluster_args = list(too_long = 1:10)
    )
  )
})





test_that("extract_pc_matrix works as expected", {
  pc_mat_sce <- extract_pc_matrix(sce)
  expect_identical(
    pc_mat_sce,
    test_mat
  )

  pc_mat_srt <- extract_pc_matrix(srat)
  expect_identical(
    pc_mat_srt,
    test_mat
  )
})
