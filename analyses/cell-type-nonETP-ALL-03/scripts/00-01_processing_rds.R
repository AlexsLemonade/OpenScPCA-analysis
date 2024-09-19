#sudo apt install libglpk40
#sudo apt install libcurl4-openssl-dev 
#running above two commands before install.packages("Seurat")
library(Seurat)
library(SingleCellExperiment) #BiocManager::install("SingleCellExperiment")
library(dplyr)
library(cowplot)
library(ggplot2)
library(future)
library(purrr)
######## step01 #######################################
library(reticulate)
#sudo apt-get install libhdf5-dev
#remotes::install_github("mojaveazure/seurat-disk")
library(SeuratDisk)
use_condaenv("openscpca")
os <- import("os")
anndata <- import("anndata") #conda install conda-forge::anndata
samalg <- import("samalg") #pip install sam-algorithm

out_loc <- "~/OpenScPCA-analysis/analyses/cell-type-nonETP-ALL-03/"
data_loc <- "~/OpenScPCA-analysis/data/2024-08-22/SCPCP000003/"
doublet_loc <- "~/OpenScPCA-analysis/data/2024-08-22/results/doublet-detection/SCPCP000003/"

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
sampleID <- metadata$scpca_sample_id[which(metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia")]
libraryID <- metadata$scpca_library_id[which(metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia")]

options(Seurat.object.assay.version = 'v3')
for (i in 1:length(sampleID)){
  sce <- readRDS(file.path(data_loc,sampleID[i],paste0(libraryID[i],"_processed.rds")))
  seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
  seu@assays[["RNA"]]@data <- logcounts(sce)
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  #detecting doublets
  doublets.data <- read.table(file = file.path(doublet_loc,sampleID[i],paste0(libraryID[i],"_processed_scdblfinder.tsv")),sep = "\t", header = T)
  idx <- match(colnames(seu), doublets.data$barcodes)
  seu$doublet_class <- doublets.data$class[idx]
  seu <- subset(seu, subset = percent.mt < 25)
  
  #step 01 feature selection/dimensionality using SAM
  SaveH5Seurat(seu, filename = paste0(sampleID[i],".h5Seurat"), overwrite = T)
  Convert(paste0(sampleID[i],".h5Seurat"), dest = "h5ad", overwrite = T)
  
  sam = samalg$SAM()
  sam$load_data(paste0(sampleID[i],".h5ad"))
  sam$preprocess_data()
  sam$run()
  #leiden clustering is the default algorithm in SAM
  sam$clustering(method = "leiden") #conda install conda-forge::leidenalg
  sam$save_anndata(paste0("sam.",sampleID[i],".h5ad"))
  Convert(paste0("sam.",sampleID[i],".h5ad"), dest = "h5seurat", overwrite = T)
  
  final.obj <- LoadH5Seurat(paste0("sam.",sampleID[i],".h5seurat"), meta.data = F, misc = F)
  final.obj <- AddMetaData(final.obj, seu@meta.data)
  # create a new assay to store ADT information
  final.obj[["ADT"]] <- CreateAssayObject(counts = counts(altExp(sce)))
  final.obj@assays[["ADT"]]@data <- logcounts(altExp(sce))
  tmp.obj <- anndata$read_h5ad(paste0("sam.",sampleID[i],".h5ad"))
  final.obj$leiden_clusters <- tmp.obj$obs[["leiden_clusters"]]
  saveRDS(final.obj, file.path(out_loc,"scratch",paste0(sampleID[i],".rds")))
  os$remove(paste0("sam.",sampleID[i],".h5ad"))
  os$remove(paste0("sam.",sampleID[i],".h5seurat"))
  os$remove(paste0(sampleID[i],".h5ad"))
  os$remove(paste0(sampleID[i],".h5Seurat"))
  
  DimPlot(final.obj, reduction = "umap", group.by = "leiden_clusters", label = T) + 
    ggtitle(paste0(sampleID[i],": leiden_clusters"))
  ggsave(file.path(out_loc,"plots/00-01_processing_rds/",paste0(sampleID[i],".pdf")))
}

