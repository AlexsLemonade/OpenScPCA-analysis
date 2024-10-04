#!/usr/bin/env Rscript

#This script combines UMAP plots of all samples into one multi-panel plot

library(Seurat)
library(ggplot2)

multiplot <- function(meta.variable, n_row=3){
  if (meta.variable %in% c("sctype_classification","lowConfidence_annot")){
    clrs <- ct_color
  }else{
    clrs <- NULL
  }
  plot.list <- list()
  for (i in 1:length(libraryID)){
    tryCatch({
      plot.list[[i]] <- DimPlot(seu.list[[i]], reduction = "Xumap_", group.by = meta.variable, 
                                label = T, cols = clrs, repel = T) +
        ggtitle(libraryID[i]) + NoLegend() + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
    }, error=function(e){})
  }
  cowplot::plot_grid(plotlist = plot.list, nrow = n_row)
  ggsave(file.path(out_loc,"plots",paste0("multipanels_",meta.variable,".png")), width = 20, height = 12, bg = "white")
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id

ct_color <- c("darkorchid","skyblue2","dodgerblue2","gold","beige","sienna1","green4","navy",
              "chocolate4","red","darkred","#6A3D9A","maroon","yellow4","grey35","black","lightpink","grey80")
names(ct_color) <- c("B","CD4 T","CD8 T","DC","HSPC","Mono","NK","Other T","Macrophage",
                     "Early Eryth","Late Eryth","Plasma","Platelet","Stromal","Blast","Cancer","Pre Eryth","Unknown")

seu.list <- list()
for (i in 1:length(libraryID)){
  seu.list[[i]] <- readRDS(file.path(out_loc,"results/rds",paste0(libraryID[i],".rds")))
}
voi <- c("leiden_clusters","sctype_classification","lowConfidence_annot","copykat.pred")

purrr::walk(voi, multiplot)
