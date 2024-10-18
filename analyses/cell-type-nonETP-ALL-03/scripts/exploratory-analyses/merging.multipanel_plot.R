#!/usr/bin/env Rscript

#This script combines multiple UMAP plots from an individual library into one multi-panel plot

library(Seurat)
library(ggplot2)

multiplot <- function(annot.obj, library.id, colors, n.row = 2, 
                      variables.to.plot = c("seurat_clusters","sctype_classification","libraryID","copykat.pred")){
  plot.list <- list()
  for (plot.type in variables.to.plot){
    if (plot.type == "sctype_classification"){
      clrs <- colors[[1]]
    } else if (plot.type == "libraryID"){
      clrs <- colors[[2]]
    } else{
      clrs <- NULL
    }
    
    plot.list[[plot.type]] <- DimPlot(annot.obj, reduction = "umap", group.by = plot.type, 
                                      label = T, cols = clrs, repel = T) + NoLegend()
  }
  cowplot::plot_grid(plotlist = plot.list, nrow = n.row) + 
    cowplot::draw_figure_label(library.id, position = "top", size = 18, fontface = "bold")
  ggsave(file.path(out_loc,"plots",paste0("multipanels_",library.id,".png")), width = 12, height = 12, bg = "white", dpi = 150)
}

multi_splitPlot <- function(annot.obj, library.id, colors){
  p1 <- DimPlot(annot.obj, reduction = "umap", group.by = "leiden_clusters",  
                split.by = "libraryID", cols = colors[[2]])
  p2 <- DimPlot(annot.obj, reduction = "umap", group.by = "sctype_classification",  
                split.by = "libraryID", cols = colors[[1]])
  cowplot::plot_grid(plotlist = list(p1,p2), nrow = 2)
  ggsave(file.path(out_loc,"plots",paste0(library.id,"_splitPlot.png")), width = 12, height = 12, bg = "white", dpi = 150)
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
lib_color <- c(ggsci::pal_jco("default",alpha = 0.4)(10),"purple")
names(lib_color) <- sort(libraryID)
c30 <- c("dodgerblue2", "#E31A1C","green4","#6A3D9A", "#FF7F00","orchid1", "gold1","skyblue2", "#FB9A99", # lt pink
         "palegreen2","#CAB2D6", # lt purple
         "darkorange4", "brown","#FDBF6F", # lt orange
         "gray70", "khaki2","maroon", "black", "deeppink1", "blue1", "steelblue4","darkturquoise", 
         "green1", "yellow4", "yellow3",
         "beige","cyan","darkgreen","navy","darkorchid")

merged.libraryID <- c("SCPCL000077-SCPCL000703","SCPCL000082-SCPCL000703","SCPCL000704-SCPCL000703",
                      "SCPCL000706-SCPCL000703","SCPCL000710-SCPCL000703")                      
seu.list <- list()
for (lib_iter in 1:length(merged.libraryID)){
  seu.list[[lib_iter]] <- readRDS(file.path(out_loc,"results/rds",paste0(merged.libraryID[lib_iter],".rds")))
  names(seu.list)[lib_iter] <- merged.libraryID[lib_iter]
}

purrr::walk2(seu.list, names(seu.list), 
             ~ multiplot(annot.obj = .x, library.id = .y, colors = list(ct_color, lib_color)))
purrr::walk2(seu.list, names(seu.list),
             ~ multi_splitPlot(annot.obj = .x, library.id = .y, colors = list(ct_color, c30)))
