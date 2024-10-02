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
  sam$save_anndata(paste0("sam.",ind.sample,".h5ad"))
  
  final.obj <- schard::h5ad2seurat(paste0("sam.",ind.sample,".h5ad"))
  final.obj <- AddMetaData(final.obj, as.data.frame(colData(sce)))
  final.obj[["RNA"]]@counts <- counts(sce)
  final.obj[["RNA"]] <- AddMetaData(final.obj[["RNA"]], as.data.frame(rowData(sce)))
  final.obj$nCount_RNA <- final.obj$sum
  final.obj$nFeature_RNA <- final.obj$detected
  # create a new assay to store ADT information
  final.obj[["ADT"]] <- CreateAssayObject(counts = counts(altExp(sce, "adt")))
  final.obj[["ADT"]]@data <- logcounts(altExp(sce, "adt"))
  final.obj[["ADT"]] <- AddMetaData(final.obj[["ADT"]], as.data.frame(rowData(altExp(sce, "adt"))))
  
  saveRDS(final.obj, file.path(out_loc,"results/rds",paste0(ind.sample,".rds")))
  unlink(paste0("sam.",ind.sample,".h5ad"))

  DimPlot(final.obj, reduction = "Xumap_", group.by = "leiden_clusters", label = T) + 
    ggtitle(paste0(ind.sample,": leiden_clusters"))
  ggsave(file.path(out_loc,"plots",paste0(ind.sample,"_leidenCluster.pdf")), width = 7, height = 7)
}

run_multisam <- function(ind.sam, libraries){
  seu.list <- list()
  rowData_RNA <- c()
  rowData_ADT <- c()
  for (i in 1:length(libraries)){
    seu.list[[i]] <- as.Seurat(readRDS(file.path(data_loc,ind.sam,paste0(libraries[i],"_processed.rds"))))
    #reading in doublets information
    doublets.data <- read.table(file = file.path(doublet_loc,ind.sam,paste0(libraries[i],"_processed_scdblfinder.tsv")),
                                sep = "\t", header = T)
    idx <- match(colnames(seu.list[[i]]), doublets.data$barcodes)
    seu.list[[i]]$doublet_class <- doublets.data$class[idx]
    seu.list[[i]]$libraryID <- libraries[i]
    if (i == 1){
      rowData_RNA <- as.data.frame(seu.list[[i]][["originalexp"]]@meta.features)
      rowData_ADT <- as.data.frame(seu.list[[i]][["adt"]]@meta.features)
    }
    else{
      rowData_RNA <- cbind(rowData_RNA, as.data.frame(seu.list[[i]][["originalexp"]]@meta.features[,3:4]))
      rowData_ADT <- cbind(rowData_ADT, as.data.frame(seu.list[[i]][["adt"]]@meta.features[,3:4]))
    }
  }
  colnames(rowData_RNA)[3:ncol(rowData_RNA)] <- paste0(rep(libraries, each = 2),"-",colnames(rowData_RNA)[3:ncol(rowData_RNA)]) 
  colnames(rowData_ADT)[3:ncol(rowData_ADT)] <- paste0(rep(libraries, each = 2),"-",colnames(rowData_ADT)[3:ncol(rowData_ADT)]) 
  sample.combined <- merge(seu.list[[1]], y = seu.list[2:length(libraries)], add.cell.ids = libraries, project = ind.sam)
  sample.combined$nCount_originalexp <- NULL
  sample.combined$nFeature_originalexp <- NULL
  DefaultAssay(sample.combined) <- "originalexp"
  
  #step 01 feature selection/dimensionality using SAM
  sam = samalg$SAM(counts = c(r_to_py(t(sample.combined[["originalexp"]]@counts)), 
                              r_to_py(as.array(rownames(sample.combined))),
                              r_to_py(as.array(colnames(sample.combined)))))
  sam$preprocess_data()
  sam$run(distance = 'correlation')
  sam$clustering(method = "leiden") #leiden clustering is the default algorithm in SAM
  sam$save_anndata(paste0("sam.",ind.sam,".h5ad"))
  
  final.obj <- schard::h5ad2seurat(paste0("sam.",ind.sam,".h5ad"))
  final.obj <- AddMetaData(final.obj, as.data.frame(sample.combined@meta.data))
  final.obj[["RNA"]]@counts <- sample.combined[["originalexp"]]@counts
  final.obj[["RNA"]] <- AddMetaData(final.obj[["RNA"]], rowData_RNA)
  final.obj$nCount_RNA <- final.obj$sum
  final.obj$nFeature_RNA <- final.obj$detected
  # create a new assay to store ADT information
  final.obj[["ADT"]] <- CreateAssayObject(counts = sample.combined[["adt"]]@counts)
  final.obj[["ADT"]]@data <- sample.combined[["adt"]]@data
  final.obj[["ADT"]] <- AddMetaData(final.obj[["ADT"]], rowData_ADT)
  
  saveRDS(final.obj, file.path(out_loc,"results/rds",paste0(ind.sam,".rds")))
  unlink(paste0("sam.",ind.sam,".h5ad"))
  
  DimPlot(final.obj, reduction = "Xumap_", group.by = "leiden_clusters", split.by = "libraryID", label = T) + 
    ggtitle(paste0(ind.sam,": leiden_clusters"))
  ggsave(file.path(out_loc,"plots",paste0(ind.sam,"_leidenCluster.pdf")), width = 10, height = 5)
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-ETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)
doublet_loc <- file.path(project_root, "data/current/results/doublet-detection",projectID)

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
sampleID <- metadata$scpca_sample_id
libraryID <- metadata$scpca_library_id

## removing SCPCS000066 from the list of sampleID, because there are two libraries for this sample
special.case <- c("SCPCS000066")
special.library <- libraryID[which(sampleID==special.case)]
sampleID <- setdiff(sampleID, special.case)
libraryID <- setdiff(libraryID, special.library)

purrr::walk2(sampleID, libraryID, run_sam)

run_multisam(special.case, special.library)


