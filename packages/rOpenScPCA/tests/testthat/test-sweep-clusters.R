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

  expect_length(sweep_list, 4)

  # check colnames one at a time
  sweep_list |>
    purrr::walk(
      \(df) {
        expect_setequal(
          colnames(df),
          c("cell_id", "cluster", "algorithm", "weighting", "nn", "resolution")
        )

        expect_equal(df$cell_id, colnames(sce))

        expect_s3_class(
          df$cluster,
          "factor"
        )

        expect_equal(
          unique(df$algorithm),
          "louvain"
        )
        expect_equal(
          unique(df$weighting),
          "jaccard"
        )
        expect_true(
          all(df$nn == 10) || all(df$nn == 15)
        )
        expect_true(
          all(df$resolution == 0.5) || all(df$resolution == 1)
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

        expect_equal(
          unique(df$algorithm),
          "leiden"
        )

        expect_equal(
          unique(df$objective_function),
          "modularity"
        )

        expect_true(
          all(unique(df$resolution) == 0.5) || all(unique(df$resolution) == 1)
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

  # first count algorithms
  alg_count_df <- sweep_list |>
    dplyr::bind_rows(.id = "id") |>
    dplyr::group_by(id) |>
    dplyr::summarize(alg = unique(algorithm))

  expect_equal(
    nrow(alg_count_df),
    6
  )

  expect_equal(
    sum(alg_count_df$alg == "louvain"),
    4
  )

  expect_equal(
    sum(alg_count_df$alg == "walktrap"),
    2
  )

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
