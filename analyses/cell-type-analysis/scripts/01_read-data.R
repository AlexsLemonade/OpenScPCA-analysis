library(SingleCellExperiment)

sce <- readRDS("../../data/current/SCPCL999990_processed.rds")

library(dplyr)

coldata_df <- colData(sce) |>
  as.data.frame() |>
  dplyr::select(barcodes,
                cluster,
                singler_celltype_annotation,
                cellassign_celltype_annotation)

readr::write_tsv("results/sce_annotations.tsv")
