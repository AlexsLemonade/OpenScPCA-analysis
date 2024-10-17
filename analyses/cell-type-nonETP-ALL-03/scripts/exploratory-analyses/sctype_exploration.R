#!/usr/bin/env Rscript

## This script investigates how solid are the B cells annotation, by checking the ScType score and the purity of that cluster

library(ggridges)
library(ggplot2)
library(dplyr)

Bcell_check <- function(ind.lib, methods = c("ScType","SingleR","CellAssign"),
                        variables.to.plot = c("sctype_classification","singler_celltype_annotation","cellassign_celltype_annotation")){
  seu <- readRDS(file.path(out_loc,"results/rds",paste0(ind.lib,".rds")))
  plot.list <- list()
  for (var_iter in 1:length(variables.to.plot)){
    plot.type <- variables.to.plot[var_iter]
    if (plot.type == "sctype_classification"){
      Bcell.names <- colnames(seu)[which(seu$sctype_classification=="B")]
    } else{
      celltype <- unique(seu@meta.data[[plot.type]])
      seu_df <- data.frame(FetchData(seu, vars = plot.type))|> tibble::rownames_to_column(var = "cell_id")
      seu_df <- seu_df[which(seu_df[[plot.type]] %in% celltype[grep("B cell",celltype)]),]
      Bcell.names <- seu_df$cell_id
    }
    
    if (length(Bcell.names) == 0){next}
    sctype.score <- read.table(file.path(out_loc,"results",paste0(ind.lib,"_sctype_scores.txt")), 
                               sep = "\t", header = T)
    df <- sctype.score[match(Bcell.names, rownames(sctype.score)),] %>% 
      tidyr::pivot_longer(cols = colnames(sctype.score), names_to = "celltype", values_to = "ScType.score")
    df$celltype <- gsub("\\."," ", df$celltype)
    p1 <- ggplot(df, aes(x = ScType.score, y = forcats::fct_reorder(celltype,ScType.score), fill = celltype)) + 
      geom_density_ridges() + theme_ridges() + 
      theme(legend.position = "none", axis.title.x = element_text(hjust=0.5), axis.title.y = element_text(vjust=0.5)) + 
      scale_fill_manual(values = ct_color) + xlab ("ScType score") + 
      ylab(expr(bold(!!methods[var_iter])*~"("*!!length(Bcell.names)*")"))

    eval.clus <- read.table(file.path(out_loc,"results/evalClus/",paste0(ind.lib,"_sil-purity_perCell.txt")), 
                            sep = "\t", header = T)
    plot.df <- data.frame(cell_id=Bcell.names,
                          purity=eval.clus$purity[match(Bcell.names, eval.clus$cell_id)],
                          sctypeScore=sctype.score$B[match(Bcell.names, rownames(sctype.score))],
                          leidenCluster=as.factor(seu$leiden_clusters[match(Bcell.names,colnames(seu))]))
    p2 <- ggplot(plot.df, aes(x = purity, y = sctypeScore, color = leidenCluster)) + 
      geom_point(size = 0.5) + theme_classic() + ylab("B cell ScType score") + xlim(0,1)
    
    plot.list <- c(plot.list, list(p1, p2))
  }
  if (length(plot.list) == 0){return()}
  cowplot::plot_grid(plotlist = plot.list, nrow = 3) + 
    patchwork::plot_annotation(title = paste0(ind.lib,": B cells identified in different methods")) &  
    theme(plot.title = element_text(hjust = 0.5, face="bold")) 
  ggsave(file.path(out_loc,"plots/sctype_exploration",paste0(ind.lib,"_Bcells.png")), 
         width = 10, height = 15, bg = "white", dpi = 150)
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)
dir.create(file.path(out_loc, "plots/sctype_exploration"), showWarnings = FALSE)

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id

ct_color <- c("darkorchid","skyblue2","dodgerblue2","gold","beige","sienna1","green4","navy",
              "chocolate4","red","darkred","#6A3D9A","maroon","yellow4","grey35","black","lightpink","grey80")
names(ct_color) <- c("B","CD4 T","CD8 T","DC","HSPC","Mono","NK","Other T","Macrophage",
                     "Early Eryth","Late Eryth","Plasma","Platelet","Stromal","Blast","Cancer","Pre Eryth","Unknown")

purrr::walk(libraryID, Bcell_check)
