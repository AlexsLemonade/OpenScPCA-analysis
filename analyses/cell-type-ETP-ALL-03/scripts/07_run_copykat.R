#!/usr/bin/env Rscript

library(Seurat)
library(copykat)
library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--testing"),
    action = "store_true",
    default = FALSE,
    help = "Whether we are running in test mode; KS.cut will be set to 0.025 if TRUE."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))
if (opts$testing) {
  ks.cut_param <- 0.025 # smaller value for test data
  copykat::copykat_options(KS.cut = 0.025)
} else {
  ks.cut_param <- 0.1 # default
}

run_copykat <- function(ind.lib, ks_cut = 0.1) {
  seu <- readRDS(file.path(out_loc, "results/rds", paste0(ind.lib, ".rds")))
  annot.file <- file.path(out_loc, "results", paste0(ind.lib, "_newB-normal-annotation.txt"))
  if (file.exists(annot.file)) { # the sample has new B cells annotated
    annot.df <- read.table(annot.file,
      header = F, row.names = 1, sep = "\t", stringsAsFactors = FALSE,
      colClasses = c("character", "character")
    )
    norm.cells <- rownames(annot.df)[which(annot.df$V2 == "new B")]
    n_cores <- parallel::detectCores() - 1
    copykat.test <- copykat(
      rawmat = seu@assays[["RNA"]]@counts, id.type = "Ensemble",
      ngene.chr = 5, win.size = 25, KS.cut = ks_cut, sam.name = ind.lib,
      distance = "euclidean", norm.cell.names = norm.cells,
      output.seg = "FALSE", plot.genes = "TRUE", genome = "hg20", n.cores = n_cores
    )
    idx <- match(colnames(seu), copykat.test$prediction$cell.names)
    seu$newB.copykat.pred <- copykat.test$prediction$copykat.pred[idx]
    saveRDS(seu, file = file.path(out_loc, "results/rds", paste0(ind.lib, ".rds")))
  }
}

project_root <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-ETP-ALL-03")
data_loc <- file.path(project_root, "data/current", projectID)
setwd(file.path(out_loc, "results/copykat_output"))

metadata <- read.table(file.path(data_loc, "single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
  metadata$diagnosis == "Early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id

purrr::walk(libraryID, run_copykat, ks_cut = ks.cut_param)
