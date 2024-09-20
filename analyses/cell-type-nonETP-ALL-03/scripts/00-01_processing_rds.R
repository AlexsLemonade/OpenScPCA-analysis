library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(cowplot)
library(ggplot2)
library(future)
library(purrr)
######## step01 #######################################
library(reticulate)
use_condaenv("openscpca-cell-type-nonETP-ALL-03")
os <- import("os")
anndata <- import("anndata")
samalg <- import("samalg")

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current/SCPCP000003")
doublet_loc <- file.path(project_root, "data/current/results/doublet-detection/SCPCP000003")

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
sampleID <- metadata$scpca_sample_id[which(metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia")]
libraryID <- metadata$scpca_library_id[which(metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia")]

options(Seurat.object.assay.version = 'v3')
for (i in 1:length(sampleID)){
  sce <- readRDS(file.path(data_loc,sampleID[i],paste0(libraryID[i],"_processed.rds")))
  seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
  seu@assays[["RNA"]]@data <- logcounts(sce)
  # create a new assay to store ADT information
  seu[["ADT"]] <- CreateAssayObject(counts = counts(altExp(sce)))
  seu[["ADT"]]@data <- logcounts(altExp(sce))
  mito.features <- rownames(sce)[which(grepl("^MT-", rowData(sce)$gene_symbol))]
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito.features)
  
  #detecting doublets
  doublets.data <- read.table(file = file.path(doublet_loc,sampleID[i],paste0(libraryID[i],"_processed_scdblfinder.tsv")),sep = "\t", header = T)
  idx <- match(colnames(seu), doublets.data$barcodes)
  seu$doublet_class <- doublets.data$class[idx]
  seu <- subset(seu, subset = percent.mt < 25)
  
  #step 01 feature selection/dimensionality using SAM
  sam = samalg$SAM(counts = c(r_to_py(t(seu@assays[["RNA"]]@counts)), r_to_py(as.array(rownames(seu))), 
                             r_to_py(as.array(colnames(seu)))))
  sam$preprocess_data()
  sam$run()
  sam$clustering(method = "leiden") #leiden clustering is the default algorithm in SAM
  sam$save_anndata(paste0("sam.",sampleID[i],".h5ad"))
  
  final.obj <- schard::h5ad2seurat(paste0("sam.",sampleID[i],".h5ad"))
  final.obj <- AddMetaData(final.obj, seu@meta.data)
  final.obj@assays[["RNA"]]@counts <- seu@assays[["RNA"]]@counts
  # create a new assay to store ADT information
  final.obj[["ADT"]] <- CreateAssayObject(counts = seu@assays[["ADT"]]@counts)
  final.obj@assays[["ADT"]]@data <- seu@assays[["ADT"]]@data
  saveRDS(final.obj, file.path(out_loc,"scratch",paste0(sampleID[i],".rds")))
  unlink(paste0("sam.",sampleID[i],".h5ad"))

  DimPlot(final.obj, reduction = "Xumap_", group.by = "leiden_clusters", label = T) + 
    ggtitle(paste0(sampleID[i],": leiden_clusters"))
  ggsave(file.path(out_loc,"plots/00-01_processing_rds",paste0(sampleID[i],".pdf")))
}

