library(dplyr)
library(Seurat)
library(ggpubr)
library(zellkonverter)
library(SingleCellExperiment)

prepare_fetal_atlas <- function(in_fetal_atlas = in_fetal_atlas,
                                out_fetal_atlas = out_fetal_atlas,
                                use_exist = T){
  path_h5ad <- file.path(in_fetal_atlas)
  path_out <- file.path(out_fetal_atlas)
  # do not re-prepare reference atlas
  if (file.exists(path_out) & use_exist) {
    print("fetal atlas already exist")
    return(0)
  }
  
  # download kidney fetal atlas, do not re-download
  # if (!file.exists(path_h5ad)) {
  #   download.file('https://cellgeni.cog.sanger.ac.uk/kidneycellatlas/Fetal_full_v3.h5ad', 
  #                 destfile = path_h5ad, 
  #                 method = "wget")
  # }
  
  # read h5ad into SCE here, to simplify renaming features to ensembl IDs
  sce <- zellkonverter::readH5AD(path_h5ad)
  rownames(sce) <- SingleCellExperiment::rowData(sce)$ID
  seurat_obj <- SeuratObject::CreateSeuratObject(counts = SingleCellExperiment::counts(sce),
                                                 assay = "RNA",
                                                 project = "kidneyatlas")
  # convert colData and rowData to data.frame for use in the Seurat object
  cell_metadata <- as.data.frame(SingleCellExperiment::colData(sce))
  row_metadata <- as.data.frame(SingleCellExperiment::rowData(sce))
  # add cell metadata (colData) from SingleCellExperiment to Seurat
  seurat_obj@meta.data <- cell_metadata
  # add row metadata (rowData) from SingleCellExperiment to Seurat
  seurat_obj[["RNA"]]@meta.data <- row_metadata
  # add metadata from SingleCellExperiment to Seurat
  seurat_obj@misc <- S4Vectors::metadata(sce)
  # log transform counts
  seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  ndims <- 50
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = ndims)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:ndims)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.8, algorithm = 1)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:ndims)
  # Seurat::DimPlot(seurat_obj, reduction = "umap", label = T, group.by = "compartment")
  SeuratObject::SaveSeuratRds(seurat_obj, file = path_out)
  
}