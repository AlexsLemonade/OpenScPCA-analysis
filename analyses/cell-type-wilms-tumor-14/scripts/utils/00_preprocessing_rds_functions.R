library(dplyr)
library(Seurat)
library(ggpubr)
pre_seuratobj <- function(obj, nfeatures = 500, run_harmony = TRUE, reduction = "harmony", ndims = 50,
                          skip_logNorm = TRUE){
  ######## Normalize, scale, feature selection
  if (skip_logNorm){
    # double check if data layer is present or not
    stopifnot("data layer must be present if skip_logNorm selected" = 
      "data" %in% SeuratObject::Layers(obj[["RNA"]])
    )
  }else{
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize")
  }
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

## create and process a seurat object from provided sce object
process_sce <- function( sample, library,
                         path_proj,
                         db_proj,
                         path_anal,
                         get_logcounts = TRUE
){
  
  ########### per sample rds preprocessing ##########
  
  #sample <- "SCPCS000517"; library <- "SCPCL000849"
  rds <- readRDS( file.path(path_proj, sample, paste0(library,"_processed.rds")) )
  db <- read.table(file.path(db_proj, sample, paste0(library,"_processed_scdblfinder.tsv")), header = T)
  

  ######## annotate doublets
  #rds <- rds[,which(db$class == "singlet")]
  rds$doublet_class <- db$class
  
  ######## create seurat object from the SCE counts matrix
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(rds),assay = "RNA",project = library)
  
  # add normalized data if requested
  if (get_logcounts) {
    SeuratObject::LayerData(seurat_obj[["RNA"]], "data") <- SingleCellExperiment::logcounts(rds)
  }
  # convert colData and rowData to data.frame for use in the Seurat object
  cell_metadata <- as.data.frame(SingleCellExperiment::colData(rds))
  row_metadata <- as.data.frame(SingleCellExperiment::rowData(rds))
  # add cell metadata (colData) from SingleCellExperiment to Seurat
  seurat_obj@meta.data <- cell_metadata
  # add row metadata (rowData) from SingleCellExperiment to Seurat
  seurat_obj[["RNA"]]@meta.data <- row_metadata
  # add metadata from SingleCellExperiment to Seurat
  seurat_obj@misc <- S4Vectors::metadata(rds)
  
  seurat_obj <- pre_seuratobj(seurat_obj, nfeatures = 3000, run_harmony = F, reduction = "pca", skip_logNorm = get_logcounts)
  
  # Seurat::DimPlot(seurat_obj, reduction = "umap", label = T)
  # obj <- Seurat::RunTSNE(obj, dims = 1:ndims)
  # Seurat::DimPlot(obj, reduction = "tsne")

  dir.create(file.path(path_anal, "scratch", "00_preprocessing_rds"),showWarnings = FALSE, recursive = TRUE)
  SeuratObject::SaveSeuratRds(seurat_obj, file = file.path(path_anal,"scratch","00_preprocessing_rds",paste0(sample,".rdsSeurat")) )
}