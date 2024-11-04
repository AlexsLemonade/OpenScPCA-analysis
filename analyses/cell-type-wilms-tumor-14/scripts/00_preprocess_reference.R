#!/usr/bin/env Rscript
# parse arguments from command line

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--in_fetal_atlas"),
    type = "character",
    default = NULL,
    help = "Path to kidney atlas (fetal) in h5ad"
  ),
  make_option(
    opt_str = c("--out_fetal_atlas"),
    type = "character",
    default = NULL,
    help = "Path to converted seurat object for kidney atlas (fetal)"
  ),
  make_option(
    opt_str = c("--run_SCT"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Generate the SCT normalized atlas"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

# make sure all input files exist
stopifnot(
  "in_fetal_atlas does not exist" = file.exists(opt$in_fetal_atlas)
)

########### function for reference preparation ##########
library(dplyr)
library(Seurat)
library(ggpubr)
library(zellkonverter)
library(SingleCellExperiment)
library(glmGamPoi)

prepare_fetal_atlas <- function(in_fetal_atlas,
                                out_fetal_atlas,
                                use_exist = TRUE) {
  # do not re-prepare reference atlas
  if (file.exists(out_fetal_atlas) & use_exist) {
    print("fetal atlas already exist")
    return(invisible(out_fetal_atlas))
  }

  # download kidney fetal atlas, do not re-download
  # if (!file.exists(path_h5ad)) {
  #   download.file('https://cellgeni.cog.sanger.ac.uk/kidneycellatlas/Fetal_full_v3.h5ad',
  #                 destfile = path_h5ad,
  #                 method = "wget")
  # }

  # read h5ad into SCE here, to simplify renaming features to ensembl IDs
  sce <- zellkonverter::readH5AD(in_fetal_atlas)
  rownames(sce) <- SingleCellExperiment::rowData(sce)$ID
  seurat_obj <- SeuratObject::CreateSeuratObject(
    counts = SingleCellExperiment::counts(sce),
    assay = "RNA",
    project = "kidneyatlas"
  )
  # convert colData and rowData to data.frame for use in the Seurat object
  cell_metadata <- as.data.frame(SingleCellExperiment::colData(sce))
  row_metadata <- as.data.frame(SingleCellExperiment::rowData(sce))
  # add cell metadata (colData) from SingleCellExperiment to Seurat
  seurat_obj@meta.data <- cell_metadata
  # add row metadata (rowData) from SingleCellExperiment to Seurat
  seurat_obj[["RNA"]]@meta.data <- row_metadata
  # add metadata from SingleCellExperiment to Seurat
  seurat_obj@misc <- S4Vectors::metadata(sce)

  # normalization
  options(future.globals.maxSize = 2000 * 1024^2)
  # log transform counts using strandard seurat workflow
  seurat_obj <- Seurat::NormalizeData(seurat_obj, normalization.method = "LogNormalize")
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  # normalize with SCTransform if requested
  if (opt$run_SCT) {
    seurat_obj <- Seurat::SCTransform(seurat_obj, conserve.memory = TRUE)
  }

  ndims <- 50
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = ndims)
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:ndims)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.8, algorithm = 1)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:ndims)
  # Seurat::DimPlot(seurat_obj, reduction = "umap", label = T, group.by = "compartment")
  SeuratObject::SaveSeuratRds(seurat_obj, file = out_fetal_atlas)
}

########### Prepare reference seurat obj ##########
prepare_fetal_atlas(
  in_fetal_atlas = opt$in_fetal_atlas,
  out_fetal_atlas = opt$out_fetal_atlas,
  use_exist = TRUE
)
