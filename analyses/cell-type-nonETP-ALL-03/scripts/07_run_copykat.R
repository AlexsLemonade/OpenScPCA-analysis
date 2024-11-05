#!/usr/bin/env Rscript

library(Seurat)
library(copykat)

run_copykat <- function(ind.lib){
  seu <- readRDS(file.path(out_loc,"results/rds",paste0(ind.lib,".rds")))
  annot.file <- file.path(out_loc,"results",paste0(ind.lib,"_newB-normal-annotation.txt"))
  if (file.exists(annot.file)){  #the sample has new B cells annotated
    annot.df <- read.table(annot.file, header=F, row.names=1, sep="\t", stringsAsFactors=FALSE, 
                           colClasses = c('character', 'character'))
    norm.cells <- rownames(annot.df)[which(annot.df$V2=="new B")]
    n_cores <- parallel::detectCores() - 1
    copykat.test <- copykat(rawmat=seu@assays[["RNA"]]@counts, id.type="Ensemble",
                            ngene.chr=5, win.size=25, KS.cut=0.1, sam.name=ind.lib,
                            distance="euclidean", norm.cell.names=norm.cells,
                            output.seg="FALSE", plot.genes="TRUE", genome="hg20",n.cores=n_cores)
    idx <- match(colnames(seu), copykat.test$prediction$cell.names)
    seu$newB.copykat.pred <- copykat.test$prediction$copykat.pred[idx]
    saveRDS(seu, file = file.path(out_loc,"results/rds",paste0(ind.lib,".rds")))
  }
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)
setwd(file.path(out_loc,"results/copykat_output"))

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id

purrr::walk(libraryID, run_copykat)
