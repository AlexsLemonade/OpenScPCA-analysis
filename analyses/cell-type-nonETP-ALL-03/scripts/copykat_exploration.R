#!/usr/bin/env Rscript

#This script explores the results of CopyKat prediction with respective of cell types and blast module scores

library(Seurat)
library(ggplot2)
library(dplyr)

copykatInterpret <- function(annot.obj, library.id, ct.colors){
  tryCatch({
    exprs <- data.frame(FetchData(annot.obj, vars = c("sctype_classification","copykat.pred","lowConfidence_annot")))
    df <- exprs %>%
      dplyr::group_by(copykat.pred, sctype_classification) %>%
      dplyr::count(name = "proportion")
    p1 <- ggplot(df, aes(x = copykat.pred, y = proportion, fill =  sctype_classification)) +
            geom_bar(width = 0.5, stat = "identity", position = "fill")+scale_fill_manual(values = ct_color)
    
    df <- exprs %>%
      dplyr::group_by(copykat.pred, lowConfidence_annot) %>%
      dplyr::count(name = "proportion")
    p2 <- ggplot(df, aes(x = copykat.pred, y = proportion, fill =  lowConfidence_annot)) +
      geom_bar(width = 0.5, stat = "identity", position = "fill")+scale_fill_manual(values = ct_color)
    cowplot::plot_grid(plotlist = list(p1,p2), nrow = 1) + 
      cowplot::draw_figure_label(library.id, position = "top", size = 14, fontface = "bold")
    ggsave(file.path(out_loc,"plots/copykat_exploration",paste0(library.id,"_celltypeVScopykat.png")), width = 10, height = 5, bg = "white", dpi = 150)
    
    ### plotting blast module scores 
    Idents(annot.obj) <- factor(annot.obj$copykat.pred, levels = c("aneuploid","diploid","not.defined")) 
    VlnPlot(annot.obj, features = "Blast_Features1") + ggtitle(paste0(library.id,": Blast module score")) + NoLegend()
    ggsave(file.path(out_loc,"plots/copykat_exploration",paste0(library.id,"_blastModuleScore.png")), width = 6, height = 6, bg = "white", dpi = 150)
  }, error=function(e){})
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)
dir.create(file.path(out_loc, "plots/copykat_exploration"), showWarnings = FALSE)

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
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

purrr::walk2(seu.list, names(seu.list), ~ copykatInterpret (annot.obj = .x, library.id = .y, ct.colors = ct_color))
