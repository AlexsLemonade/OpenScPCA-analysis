#!/usr/bin/env Rscript

#This script uses SAM algorithm to perform soft feature selection for better separation of homogenous populations.
#So we use SAM to preprocess the data, perform feature selection/dimensionality reduction and clustering.
#The outputs are an intermediate rds file and umap plot showing leiden clustering results.

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
######## step01 #######################################
library(reticulate)
use_condaenv("openscpca-cell-type-ETP-ALL-03")
samalg <- import("samalg") #https://github.com/atarashansky/self-assembling-manifold/tree/master

run_sam <- function(ind.sample, ind.library){
  sce <- readRDS(file.path(data_loc,ind.sample,paste0(ind.library,"_processed.rds")))
  
  #reading in doublets information
  doublets.data <- read.table(file = file.path(doublet_loc,ind.sample,paste0(ind.library,"_processed_scdblfinder.tsv")),
                              sep = "\t", header = T)
  idx <- match(colnames(sce), doublets.data$barcodes)
  sce$doublet_class <- doublets.data$class[idx]
  sce <- sce[, which(sce$subsets_mito_percent < 25)]
  
  #step 01 feature selection/dimensionality using SAM
  sam = samalg$SAM(counts = c(r_to_py(t(counts(sce))), 
                              r_to_py(as.array(rownames(sce))),
                              r_to_py(as.array(colnames(sce)))))
  sam$preprocess_data()
  sam$run(distance = 'correlation')
  sam$clustering(method = "leiden") #leiden clustering is the default algorithm in SAM
  sam$save_anndata(paste0("sam.",ind.library,".h5ad"))
  
  final.obj <- schard::h5ad2seurat(paste0("sam.",ind.library,".h5ad"))
  final.obj <- AddMetaData(final.obj, as.data.frame(colData(sce)))
  final.obj[["RNA"]]@counts <- counts(sce)
  final.obj[["RNA"]] <- AddMetaData(final.obj[["RNA"]], as.data.frame(rowData(sce)))
  final.obj$nCount_RNA <- final.obj$sum
  final.obj$nFeature_RNA <- final.obj$detected
  # create a new assay to store ADT information
  final.obj[["ADT"]] <- CreateAssayObject(counts = counts(altExp(sce, "adt")))
  final.obj[["ADT"]]@data <- logcounts(altExp(sce, "adt"))
  final.obj[["ADT"]] <- AddMetaData(final.obj[["ADT"]], as.data.frame(rowData(altExp(sce, "adt"))))
  
  saveRDS(final.obj, file.path(out_loc,"results/rds",paste0(ind.library,".rds")))
  unlink(paste0("sam.",ind.library,".h5ad"))

  DimPlot(final.obj, reduction = "Xumap_", group.by = "leiden_clusters", label = T) + 
    ggtitle(paste0(ind.library,": leiden_clusters"))
  ggsave(file.path(out_loc,"plots",paste0(ind.library,"_leidenCluster.png")), width = 7, height = 7)
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-ETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)
doublet_loc <- file.path(project_root, "data/current/results/doublet-detection",projectID)
dir.create(file.path(out_loc, "results/rds"), showWarnings = FALSE)

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
sampleID <- metadata$scpca_sample_id
libraryID <- metadata$scpca_library_id

purrr::walk2(sampleID, libraryID, run_sam)
