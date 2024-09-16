library(dplyr)
library(Seurat)
library(ggpubr)
pre_seuratobj <- function(obj, nfeatures = 500, run_harmony = TRUE, reduction = "harmony", ndims = 50){
  ######## Normalize, scale, feature selection
  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  obj <- Seurat::ScaleData(obj, features = Seurat::VariableFeatures(object = obj))
  # obj <- Seurat::SCTransform(obj)
  obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(object = obj))
  # Seurat::ElbowPlot(obj, ndims = 50)
  
  ######## batch effect correction
  if (run_harmony){
    obj <- harmony::RunHarmony(obj, group.by.vars = "library_id")
  }
  
  ######## Clustering and dimentional reduction
  obj <- Seurat::FindNeighbors(obj, dims = 1:ndims, reduction = reduction)
  obj <- Seurat::FindClusters(obj, resolution = 0.8, algorithm = 1)
  obj <- Seurat::RunUMAP(obj, dims = 1:ndims, reduction = reduction)
  
  return(obj)
}

## remove duplicated genes in the list for dot plot
clean_gslist <- function(inlist, obj){
  dupgenes <- unique(unlist(inlist)[duplicated(unlist(inlist))])
  availgenes <- rownames(obj)
  for (i in 1:length(inlist)) {
    element <- inlist[[i]][inlist[[i]] %in% availgenes]
    element <- element[!element %in% dupgenes]
    inlist[[i]] <- element
  }
  return(inlist)
}

## create and process a seurat object from provided sce object
process_sce <- function( sample, library,
                         path_proj = "/home/lightsail-user/git/OpenScPCA-analysis/data/current/SCPCP000014",
                         db_proj = "/home/lightsail-user/git/OpenScPCA-analysis/data/current/results/doublet-detection/SCPCP000014",
                         path_anal = "/home/lightsail-user/git/OpenScPCA-analysis/analyses/cell-type-wilms-tumor-14/"
){
  
  ########### per sample rds preprocessing ##########
  
  #sample <- "SCPCS000517"; library <- "SCPCL000849"
  rds <- readRDS(paste0(path_proj,"/",sample,"/",library,"_processed.rds"))
  # out_rds <- paste0(path_proj,"/",sample,"/",library,"_processed_genesym.rds")
  db <- read.table(paste0(db_proj,"/",sample,"/",library,"_processed_scdblfinder.tsv"), header = T)
  
  ######## convert ensembl IDs to symbols and dedup
  rownames(rds) <- rds@rowRanges@elementMetadata@listData[["gene_symbol"]]
  rds <- singleCellTK::dedupRowNames(rds, as.rowData = F, return.list = F)
  rds <- rds[!is.na(rownames(rds)),]
  rds <- rds[!duplicated(rownames(rds)),]
  # readr::write_rds(rds, out_rds, compress = "bz2")
  
  ######## remove doublets
  rds <- rds[,which(db$class == "singlet")]
  
  ######## create seurat object from the SCE counts matrix
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(rds),
                                                 assay = "RNA",
                                                 project = library)
  # convert colData and rowData to data.frame for use in the Seurat object
  cell_metadata <- as.data.frame(SingleCellExperiment::colData(rds))
  row_metadata <- as.data.frame(SingleCellExperiment::rowData(rds))
  # add cell metadata (colData) from SingleCellExperiment to Seurat
  seurat_obj@meta.data <- cell_metadata
  # add row metadata (rowData) from SingleCellExperiment to Seurat
  seurat_obj[["RNA"]]@meta.data <- row_metadata
  # add metadata from SingleCellExperiment to Seurat
  seurat_obj@misc <- S4Vectors::metadata(rds)
  
  seurat_obj <- pre_seuratobj(seurat_obj, nfeatures = 3000, run_harmony = F, reduction = "pca")
  
  # Seurat::DimPlot(seurat_obj, reduction = "umap", label = T)
  # obj <- Seurat::RunTSNE(obj, dims = 1:ndims)
  # Seurat::DimPlot(obj, reduction = "tsne")

  dir.create(paste0(path_anal,"/results/00_preprocessing_rds"),showWarnings = FALSE, recursive = TRUE)
  SeuratObject::SaveSeuratRds(seurat_obj, file = paste0(path_anal,"/results/00_preprocessing_rds/",sample,".rdsSeurat"))
}