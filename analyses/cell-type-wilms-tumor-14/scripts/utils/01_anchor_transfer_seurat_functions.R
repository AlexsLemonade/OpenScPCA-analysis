library(dplyr)
library(Seurat)
library(ggpubr)
library(zellkonverter)
library(SingleCellExperiment)

prepare_fetal_atlas <- function(scratch_out_dir, use_exist = T){
  path_h5ad <- file.path(scratch_out_dir, "Fetal_full_v3.h5ad")
  path_out <- file.path(scratch_out_dir, "kidneyatlas.rdsSeurat")
  # do not re-prepare reference atlas
  if (file.exists(path_out) & use_exist) {
    print("fetal atlas already exist")
    return(0)
  }
  
  # download kidney fetal atlas
  if (!file.exists(path_h5ad)) {
    download.file('https://cellgeni.cog.sanger.ac.uk/kidneycellatlas/Fetal_full_v3.h5ad', 
                  destfile = path_h5ad, 
                  method = "wget")
  }
  
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
  seurat_obj <- Seurat::ScaleData(seurat_obj, features = Seurat::VariableFeatures(object = seurat_obj))
  seurat_obj <- Seurat::RunPCA(seurat_obj, features = Seurat::VariableFeatures(object = seurat_obj))
  ndims <- 50
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:ndims)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.8, algorithm = 1)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:ndims)
  # Seurat::DimPlot(seurat_obj, reduction = "umap", label = T, group.by = "compartment")
  SeuratObject::SaveSeuratRds(seurat_obj, file = path_out)
  
}

run_anchorTrans <- function(path_anal, scratch_out_dir, results_out_dir,
                            ref_obj, sample,
                            unknown_cutoff = 0.5) {
  # perspective output files
  filename <- file.path(results_out_dir, paste0(sample, "_anchor.pdf"))
  # if (file.exists(filename1) & file.exists(filename2)) {
  #   print("results already exist")
  #   return(0)
  # }
  
  sample_obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(sample,".rdsSeurat")) )

  # find anchors
  anchors <- FindTransferAnchors(reference = ref_obj, query = sample_obj)
  
  # clean annotation
  ref_obj@meta.data <- ref_obj@meta.data %>%
    mutate(annot = case_when(compartment == "stroma" ~ "stroma",
                             compartment == "immune" ~ "immune",
                             compartment == "endothelium" ~ "endothelium",
                             TRUE ~ celltype))
  ref_obj@meta.data$annot <- factor(ref_obj@meta.data$annot)
  # transfer labels
  predictions <- TransferData(
    anchorset = anchors,
    refdata = ref_obj$annot
  )
  predictions <- mutate(predictions, predicted.id = case_when(prediction.score.max < unknown_cutoff ~ "Unknown",
                                                              TRUE ~ predicted.id))
  sample_obj <- AddMetaData(object = sample_obj, metadata = predictions)
  
  # plot
  color <- Polychrome::glasbey.colors( length( unique(predictions$predicted.id) ) +1  )
  color <- color[-1]
  names(color) <- unique(predictions$predicted.id)
  color['Unknown'] <- "gray90"
  
  p1 <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", label = F, cols = color) +
    ggtitle(sample)
  p2 <- Seurat::DimPlot(sample_obj, reduction = "umap", label = T)
  df <- sample_obj@meta.data %>%
    dplyr::group_by(seurat_clusters, predicted.id) %>%
    dplyr::count(name = "sum")
  p3 <- ggplot(df, aes(x = seurat_clusters, y = sum, fill =  predicted.id)) +
    geom_bar(width = 0.5, stat = "identity", position = "fill") +
    scale_fill_manual(values = color)
  toprow <- ggpubr::ggarrange(p1, p2, widths = c(1.5,1))
  p <- ggpubr::ggarrange(toprow, p3, ncol = 1)
  p_split <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", 
                             label = F, cols = color, split.by = "predicted.id", ncol = 3, alpha = 0.1) +
    ggtitle(sample)
  
  # save plots
  multi_page <- ggpubr::ggarrange(p, p_split,
                                  nrow = 1, ncol = 1)
  ggpubr::ggexport(plotlist = multi_page,
                   filename = filename,
                   width = ifelse(length( unique(predictions$predicted.id) ) < 15,9,12), height = 6)
  
  
  
}
