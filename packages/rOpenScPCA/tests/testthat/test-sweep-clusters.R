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


test_that("sweep_clusters works as expected with default algorithm & weighting", {
  sweep_list <- sweep_clusters(
    sce,
    nn = c(10, 15),
    resolution = c(0.5, 1)
  )

  expect_length(sweep_list, 4)

  sweep_list |>
    purrr::walk(
      \(df) {
        expect_setequal(
          colnames(df),
          c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
        )

        # these tests confirm the defaults went through
        expect_equal(unique(df$algorithm), "louvain")
        expect_equal(unique(df$weighting), "jaccard")

        expect_true(
          all(df$nn == 10) || all(df$nn == 15)
        )
        expect_true(
          all(df$resolution == 0.5) || all(df$resolution == 1)
        )
      }
    )
})



test_that("sweep_clusters works as expected with matrix input", {
  sweep_list <- sweep_clusters(
    test_mat,
    nn = c(10, 15)
  )


  expect_length(sweep_list, 2)

  sweep_list |>
    purrr::walk(
      \(df) {
        expect_setequal(
          colnames(df),
          c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
        )
      }
    )
})



test_that("sweep_clusters works as expected with Seurat input", {
  sweep_list <- sweep_clusters(
    srat,
    nn = c(10, 15)
  )


  expect_length(sweep_list, 2)

  sweep_list |>
    purrr::walk(
      \(df) {
        expect_setequal(
          colnames(df),
          c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
        )
      }
    )
})


test_that("sweep_clusters works as expected with non-default algorithm", {
  sweep_list <- sweep_clusters(
    sce,
    algorithm = "leiden",
    objective_function = "modularity",
    resolution = c(0.5, 1)
  )

  sweep_list |>
    purrr::walk(
      \(df) {
        expect_setequal(
          colnames(df),
          c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution", "objective_function")
        )

        expect_equal(unique(df$algorithm), "leiden")
        expect_equal(unique(df$objective_function), "modularity")

        expect_true(
          all(df$resolution == 0.5) || all(df$resolution == 1)
        )
      }
    )
})




test_that("sweep_clusters works as expected with multiple algorithms", {
  sweep_list <- sweep_clusters(
    sce,
    algorithm = c("walktrap", "louvain"),
    # used by both
    nn = c(10, 15),
    # only used by louvain
    resolution = c(0.5, 1)
  )

  # count algorithms
  alg_counts <- sweep_list |>
    purrr::map(\(df) unique(df$algorithm)) |>
    purrr::reduce(c)
  expect_equal(length(alg_counts), 6)
  expect_equal(sum(alg_counts == "louvain"), 4)
  expect_equal(sum(alg_counts == "walktrap"), 2)



  sweep_list |>
    purrr::walk(
      \(df) {
        if (unique(df$algorithm) == "walktrap") {
          expect_setequal(
            colnames(df),
            c("cell_id", "cluster", "algorithm", "weighting", "nn")
          )
        } else if (unique(df$algorithm) == "louvain") {
          expect_setequal(
            colnames(df),
            c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
          )
        }
      }
    )
})
