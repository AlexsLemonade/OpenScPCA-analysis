#!/usr/bin/env Rscript

#This script combines multiple UMAP plots from an individual library into one multi-panel plot

library(Seurat)
library(ggplot2)

multiplot <- function(annot.obj, library.id, ct.colors, n.row = 2, 
                      variables.to.plot = c("leiden_clusters","sctype_classification","lowConfidence_annot","copykat.pred")){
  plot.list <- list()
  for (plot.type in variables.to.plot){
    if (plot.type %in% c("sctype_classification","lowConfidence_annot")){
      clrs <- ct.colors
    } else{
      clrs <- NULL
    }
    tryCatch({
      plot.list[[plot.type]] <- DimPlot(annot.obj, reduction = "Xumap_", group.by = plot.type, 
                                        label = T, cols = clrs, repel = T) + NoLegend()
    }, error=function(e){})
  }
  cowplot::plot_grid(plotlist = plot.list, nrow = n.row) + 
    cowplot::draw_figure_label(library.id, position = "top", size = 18, fontface = "bold")
  ggsave(file.path(out_loc,"plots",paste0("multipanels_",library.id,".png")), width = 12, height = 12, bg = "white", dpi = 150)
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-ETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id

ct_color <- c("darkorchid","skyblue2","dodgerblue2","gold","beige","sienna1","green4","navy",
              "chocolate4","red","darkred","#6A3D9A","maroon","yellow4","grey35","black","lightpink","grey80")
names(ct_color) <- c("B","CD4 T","CD8 T","DC","HSPC","Mono","NK","Other T","Macrophage",
                     "Early Eryth","Late Eryth","Plasma","Platelet","Stromal","Blast","Cancer","Pre Eryth","Unknown")

seu.list <- list()
for (lib_iter in 1:length(libraryID)){
  seu.list[[lib_iter]] <- readRDS(file.path(out_loc,"results/rds",paste0(libraryID[lib_iter],".rds")))
  names(seu.list)[lib_iter] <- libraryID[lib_iter]
}

purrr::walk2(seu.list, names(seu.list), ~ multiplot (annot.obj = .x, library.id = .y, ct.colors = ct_color))
