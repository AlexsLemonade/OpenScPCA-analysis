suppressPackageStartupMessages(library(SingleCellExperiment))

set.seed(2024)
sce <- splatter::simpleSimulate(nGenes = 1000, verbose = FALSE) |>
  scater::logNormCounts() |>
  scater::runPCA(ncomponents = 10)

test_mat <- reducedDim(sce, "PCA")

srat <- Seurat::CreateSeuratObject(counts = counts(sce), assay = "RNA")
srat[["pca"]] <- Seurat::CreateDimReducObject(
  embeddings = test_mat,
  key = "PC_", # underscore avoids Seurat warning that it's adding an underscore
  assay = "RNA"
)

test_that("calculate_clusters runs with a matrix, defaults", {
  cluster_df <- calculate_clusters(test_mat)

  expect_setequal(
    colnames(cluster_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
  )
  expect_equal(cluster_df$cell_id, rownames(test_mat))
  expect_s3_class(cluster_df$cluster, "factor")
  expect_equal(unique(cluster_df$algorithm), "louvain")
  expect_equal(unique(cluster_df$weighting), "jaccard")
  expect_equal(unique(cluster_df$nn), 10)
  expect_equal(unique(cluster_df$resolution), 1)
})


test_that("calculate_clusters runs with additional cluster_args", {
  cluster_df <- calculate_clusters(
    test_mat,
    algorithm = "leiden",
    cluster_args = list(n_iterations = 3)
  )

  expect_setequal(
    colnames(cluster_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution", "objective_function", "n_iterations")
  )
  expect_equal(unique(cluster_df$n_iterations), 3)
})



test_that("calculate_clusters runs when cluster_args is empty", {
  cluster_df <- calculate_clusters(
    test_mat,
    algorithm = "walktrap"
  )

  expect_setequal(
    colnames(cluster_df),
    c("cell_id", "cluster", "algorithm", "weighting", "nn")
  )
  expect_equal(unique(cluster_df$algorithm), "walktrap")
})


test_that("calculate_clusters runs with an object, defaults", {
  cluster_df_sce <- calculate_clusters(sce)
  expect_setequal(
    colnames(cluster_df_sce),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
  )
  expect_equal(cluster_df_sce$cell_id, rownames(test_mat))

  cluster_df_srat <- calculate_clusters(srat)
  expect_setequal(
    colnames(cluster_df_srat),
    c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
  )
  expect_equal(cluster_df_srat$cell_id, rownames(test_mat))
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
  # update test_mat column names to match what will have Seurat changed them to
  colnames(test_mat) <- gsub("^PC", "PC_", colnames(test_mat))
  expect_identical(pc_mat_srt, test_mat)
})

test_that("extract_pc_matrix errors as expected", {
  expect_error(
    extract_pc_matrix(sce, pc_name = "bad_name")
  )
  expect_error(
    extract_pc_matrix(srat, pc_name = "bad_name")
  )
  expect_error(
    extract_pc_matrix(test_mat)
  )
})
