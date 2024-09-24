#!/usr/bin/env Rscript

#This script uses SAM algorithm to perform soft feature selection for better separation of homogenous populations.
#So we use SAM to preprocess the data, perform feature selection/dimensionality reduction and clustering.
#The outputs are an intermediate rds file and umap plot showing leiden clustering results.

library(Seurat)
library(SingleCellExperiment)
library(ggplot2)
######## step01 #######################################
library(reticulate)
use_condaenv("openscpca-cell-type-nonETP-ALL-03")
samalg <- import("samalg") #https://github.com/atarashansky/self-assembling-manifold/tree/master

run_sam <- function(sample, library){
  sce <- readRDS(file.path(data_loc,sample,paste0(library,"_processed.rds")))
  seu <- CreateSeuratObject(counts = counts(sce), data = logcounts(sce), meta.data = as.data.frame(colData(sce)))
  # create a new assay to store ADT information
  seu[["ADT"]] <- CreateAssayObject(counts = counts(altExp(sce)))
  seu[["ADT"]]@data <- logcounts(altExp(sce))
  mito.features <- rownames(sce)[which(grepl("^MT-", rowData(sce)$gene_symbol))]
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito.features)
  
  #detecting doublets
  doublets.data <- read.table(file = file.path(doublet_loc,sample,paste0(library,"_processed_scdblfinder.tsv")),sep = "\t", header = T)
  idx <- match(colnames(seu), doublets.data$barcodes)
  seu$doublet_class <- doublets.data$class[idx]
  seu <- subset(seu, subset = percent.mt < 25)
  
  #step 01 feature selection/dimensionality using SAM
  sam = samalg$SAM(counts = c(r_to_py(t(seu@assays[["RNA"]]@layers[["counts"]])), 
                              r_to_py(as.array(rownames(seu))),
                              r_to_py(as.array(colnames(seu)))))
  sam$preprocess_data()
  sam$run()
  sam$clustering(method = "leiden") #leiden clustering is the default algorithm in SAM
  sam$save_anndata(paste0("sam.",sample,".h5ad"))
  
  final.obj <- schard::h5ad2seurat(paste0("sam.",sample,".h5ad"))
  final.obj <- AddMetaData(final.obj, seu@meta.data)
  final.obj@assays[["RNA"]]@counts <- seu@assays[["RNA"]]@layers[["counts"]]
  # create a new assay to store ADT information
  final.obj[["ADT"]] <- CreateAssayObject(counts = seu@assays[["ADT"]]@counts)
  final.obj@assays[["ADT"]]@data <- seu@assays[["ADT"]]@data
  saveRDS(final.obj, file.path(out_loc,"scratch",paste0(sample,".rds")))
  unlink(paste0("sam.",sample,".h5ad"))

  DimPlot(final.obj, reduction = "Xumap_", group.by = "leiden_clusters", label = T) + 
    ggtitle(paste0(sample,": leiden_clusters"))
  ggsave(file.path(out_loc,"plots/00-01_processing_rds",paste0(sample,".pdf")))
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)
doublet_loc <- file.path(project_root, "data/current/results/doublet-detection",projectID)

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
sampleID <- metadata$scpca_sample_id
libraryID <- metadata$scpca_library_id

purrr::walk2(sampleID, libraryID, run_sam)
