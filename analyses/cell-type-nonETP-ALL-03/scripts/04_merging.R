#!/usr/bin/env Rscript

## This script merges sample without B cells with SCPCL000710 (the best annotated sample)
## and corrects the batch effect with Harmony to see if there are any cells from the sample
## without B cells cluster with SCPCL000710, deducing their tumor states. 

library(Seurat)
library(ggplot2)
library(dplyr)

merging <- function(ind.lib, norm.lib = "SCPCL000703", num_PC = 20, res = 0.8){
  seu.list <- list()
  for (lib in c(ind.lib,norm.lib)){
    seu.list[[lib]] <- readRDS(file.path(out_loc,"results/rds",paste0(lib,".rds")))
    seu.list[[lib]]$libraryID <- lib
  }
  sample.combined <- merge(seu.list[[ind.lib]], y = seu.list[[norm.lib]], 
                           add.cell.ids = c(ind.lib, norm.lib))
  sample.combined <- NormalizeData(sample.combined, verbose = FALSE) #LogNormalize
  
  ## adding module score again because re-normalize
  gs_list <- gene_sets_prepare(db, tissue) #prepare gene sets
  for (i in 1:length(gs_list$gs_positive)){
    sample.combined <- AddModuleScore(object = sample.combined, name = paste0(gsub(" ","",names(gs_list$gs_positive[i])),"_Features"), 
                                      features = list(gs_list$gs_positive[[i]]))
  }
  final.obj <- sample.combined #avoid saving scale.data slot
  
  sample.combined <- FindVariableFeatures(sample.combined, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sample.combined <- ScaleData(sample.combined, features = rownames(sample.combined))
  sample.combined <- RunPCA(sample.combined, features = VariableFeatures(object = sample.combined))
  
  set.seed(42)
  sample.combined <- harmony::RunHarmony(sample.combined, "libraryID", plot_convergence = TRUE, nclust = 50, max_iter = 10, early_stop = T)
  sample.combined <- sample.combined %>%
    RunUMAP(reduction = "harmony", dims = 1:num_PC) %>%
    FindNeighbors(reduction = "harmony", dims = 1:num_PC) %>%
    FindClusters(resolution = res) %>%
    identity()
  final.obj$seurat_clusters <- sample.combined$seurat_clusters
  final.obj[['pca']] <- CreateDimReducObject(embeddings = sample.combined@reductions[["pca"]]@cell.embeddings, key = "PCA_", global = F, assay = "RNA")
  final.obj@reductions[["pca"]]@feature.loadings <- sample.combined@reductions[["pca"]]@feature.loadings
  final.obj@reductions[["pca"]]@stdev <- sample.combined@reductions[["pca"]]@stdev
  final.obj[['umap']] <- CreateDimReducObject(embeddings = sample.combined@reductions[["umap"]]@cell.embeddings, key = "UMAP_", global = T, assay = "RNA")
  final.obj[['harmony']] <- CreateDimReducObject(embeddings = sample.combined@reductions[["harmony"]]@cell.embeddings, key = "harmony_", global = F, assay = "RNA")
  final.obj@reductions[["harmony"]]@feature.loadings <- sample.combined@reductions[["harmony"]]@feature.loadings
  final.obj@reductions[["harmony"]]@stdev <- sample.combined@reductions[["harmony"]]@stdev
  
  saveRDS(final.obj, file.path(out_loc,"results/rds",paste0(ind.lib,"-",norm.lib,".rds")))
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)

source(file.path(out_loc, "scripts/util/gene-set-functions.R"))
db <- file.path(out_loc,"Azimuth_BM_level1.csv")
tissue <- "Immune system"  
metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id

libraryID <- libraryID[c(2,3,9,10,11)]
purrr::walk(libraryID, merging)

