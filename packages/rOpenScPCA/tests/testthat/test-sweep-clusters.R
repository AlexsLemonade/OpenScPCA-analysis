suppressPackageStartupMessages(library(SingleCellExperiment))

set.seed(2024)
sce <- splatter::simpleSimulate(nGenes = 1000, verbose = FALSE) |>
  scater::logNormCounts() |>
  scater::runPCA(ncomponents = 10)


test_that("sweep_clusters works as expected with defaults", {
  sweep_list <- sweep_clusters(
    sce,
    nn = c(10, 15),
    resolution = c(0.5, 1)
  )

  sweep_list |>
    purrr::walk(
      \(df) {
        expect_setequal(
          colnames(df),
          c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
        )
      }
    )

  check_sweep_list <- sweep_list |>
    dplyr::bind_rows()

  expect_equal(
    check_sweep_list$cell_id,
    rep(colnames(sce), 4)
  )

  expect_s3_class(
    check_sweep_list$cluster,
    "factor"
  )

  expect_equal(
    unique(check_sweep_list$algorithm),
    "louvain"
  )
  expect_equal(
    unique(check_sweep_list$weighting),
    "jaccard"
  )
  expect_setequal(
    unique(check_sweep_list$nn),
    c(10, 15)
  )
  expect_setequal(
    unique(check_sweep_list$resolution),
    c(0.5, 1)
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
      }
    )

  check_sweep_list <- sweep_list |>
    dplyr::bind_rows()


  expect_equal(
    unique(check_sweep_list$algorithm),
    "leiden"
  )

  expect_equal(
    unique(check_sweep_list$objective_function),
    "modularity"
  )

  expect_setequal(
    unique(check_sweep_list$resolution),
    c(0.5, 1)
  )
})
