library(dplyr)
library(Seurat)
library(ggpubr)
library(SingleCellExperiment)


run_anchorTrans <- function(path_anal, scratch_out_dir, results_out_dir, plots_out_dir,
                            ref_obj, library, 
                            level = "celltype", # celltype or compartment
                            k_weight = 50,
                            unknown_cutoff = 0.5, ndims = 15,
                            obj_assay = "RNA") {
  # perspective output files
  # filename <- file.path(results_out_dir, paste0(library, "_", level,".pdf"))
  filename_csv <- file.path(results_out_dir, paste0(library, "_", level,".csv"))

  
  sample_obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library,".rdsSeurat")) )
  
  # set active assay
  DefaultAssay(sample_obj) <- obj_assay
  DefaultAssay(ref_obj) <- obj_assay
  
  # set row names as gene symbol, since reference obj uses gene symbol by default
  # length(intersect(rownames(ref_obj), rownames(sample_obj))) # make sure gene symbol consistency
  # length(intersect(VariableFeatures(ref_obj), rownames(sample_obj))) # candidate genes used for transfer
  # unique <- sample_obj[["RNA"]]@meta.data$gene_ids[!duplicated(sample_obj[["RNA"]]@meta.data$gene_symbol)]
  # sample_obj <- sample_obj[unique,]
  # rownames(sample_obj) <- sample_obj[["RNA"]]@meta.data$gene_symbol
  
  # find anchors
  anchors <- FindTransferAnchors(reference = ref_obj, query = sample_obj, 
                                 dims = 1:ndims)
  nanchors <- nrow(anchors@anchors)
  # clean annotation, too few stroma, immune and endo
  # ref_obj@meta.data <- ref_obj@meta.data %>%
  #   mutate(annot = case_when(compartment == "stroma" ~ "stroma",
  #                            compartment == "immune" ~ "immune",
  #                            compartment == "endothelium" ~ "endothelium",
  #                            TRUE ~ celltype))
  # ref_obj@meta.data <- ref_obj@meta.data %>%
  #   mutate(annot = ifelse(level == "compartment", compartment, celltype))
  if (level == "compartment") {
    ref_obj@meta.data <- ref_obj@meta.data %>%
      mutate(annot = compartment)
  } else {
    ref_obj@meta.data <- ref_obj@meta.data %>%
      mutate(annot = celltype) %>%
      mutate(annot = case_when(compartment == "stroma" ~ "stroma",
                              compartment == "immune" ~ "immune",
                              compartment == "endothelium" ~ "endothelium",
                              TRUE ~ annot))
  }
  ref_obj@meta.data$annot <- factor(ref_obj@meta.data$annot)
  # transfer labels
  predictions <- TransferData(anchorset = anchors, refdata = ref_obj$annot, dims = 1:ndims, k.weight = k_weight) %>%
    mutate(predicted.id = ifelse(prediction.score.max < unknown_cutoff, "Unknown", predicted.id))
  # add prediction to metadata
  sample_obj <- AddMetaData(object = sample_obj, metadata = predictions)
  
  # save prediction table
  write.csv(predictions, file = filename_csv)
  
  # anchor transfer plots
  plot_anchorTrans(path_anal = path_anal, 
                   scratch_out_dir = scratch_out_dir, 
                   results_out_dir = results_out_dir,
                   plots_out_dir = plots_out_dir,
                   sample_obj = sample_obj,
                   library = library, 
                   level = level,
                   nanchors = nanchors)
  
}

plot_anchorTrans <- function(path_anal, scratch_out_dir, results_out_dir, plots_out_dir,
                             sample_obj,
                             library, 
                             level = "celltype",
                             nanchors = 0,
                             plot_cluster = "seurat_clusters", # by default, plot cluster generated in preprocess (SCT)
                             internal = FALSE){
  # 'internal' flag to return a core plot for future exploration
  if (internal == FALSE) {
    filename <- file.path(plots_out_dir, paste0(library, "_", level,".pdf"))
    filename_core <- file.path(plots_out_dir, paste0(library, "_", level,"_core.png"))
  }
  
  # anchor transfer plot
  pred_ids <- unique(sample_obj$predicted.id)
  color <- Polychrome::glasbey.colors( length( pred_ids ) +1  )
  color <- color[-1]
  names(color) <- pred_ids
  color['Unknown'] <- "gray90"
  
  p1 <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", 
                        label = F, cols = color, alpha = 0.1) +
    ggtitle(paste0(library))
  p2 <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = plot_cluster, label = T)
  # df <- sample_obj@meta.data %>%
  #   dplyr::group_by(seurat_clusters, predicted.id) %>%
  #   dplyr::count(name = "sum")
  p3 <- ggplot(sample_obj@meta.data, aes(x = get(plot_cluster), fill =  predicted.id)) +
    geom_bar(width = 0.5, position = "fill") +
    scale_fill_manual(values = color) +
    labs(y = "Proportion of cells")
  toprow <- ggpubr::ggarrange(p1, p2, widths = c(1.5,1))
  p <- ggpubr::ggarrange(toprow, p3, ncol = 1)
  p_split <- Seurat::DimPlot(sample_obj, reduction = "umap", group.by = "predicted.id", 
                             label = F, cols = color, split.by = "predicted.id", ncol = 3, alpha = 0.1) +
    ggtitle(library)
  p_hist <- ggplot(sample_obj@meta.data, aes(x = prediction.score.max)) +
    geom_histogram(bins = 100) +
    xlim(0,1) + 
    ggtitle(paste0("Prediction score distribution ", library, ", nanchors = ",nanchors))
  
  # for internal use
  if (internal == TRUE) return(p)
  
  # save plots
  multi_page <- ggpubr::ggarrange(p, p_split, p_hist,
                                  nrow = 1, ncol = 1)
  ggpubr::ggexport(plotlist = multi_page,
                   filename = filename,
                   width = ifelse(length( pred_ids ) < 15,9,12), 
                   height = ifelse(length( pred_ids ) < 15,6,9))
  
  # save pngs
  ggsave(p, file = filename_core,
         device = "png", dpi = 200,
         width = ifelse(length( pred_ids ) < 15,9,12),
         height = ifelse(length( pred_ids ) < 15,6,9))
}
