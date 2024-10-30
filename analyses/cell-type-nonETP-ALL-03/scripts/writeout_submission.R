#!/usr/bin/env Rscript

library(Seurat)
library(ggplot2)

writeout <- function(ind.lib, ct.colors = ct_color, project.ID = projectID, n.row = 1){
  seu <- readRDS(file.path(out_loc,"results/rds",paste0(ind.lib,".rds")))
  voi <- c('newB.copykat.pred','sctype_classification')
  changeName.voi <- c('tumor_cell_classification','cell_type_assignment')
  tryCatch({
    voi_df <- data.frame(FetchData(seu, vars = voi)) |> tibble::rownames_to_column(var = "cell_barcode")
  }, error=function(e){})
  colnames(voi_df)[2:length(voi_df)] <- changeName.voi[match(colnames(voi_df)[2:length(voi_df)],voi)]
  final.df <- data.frame(scpca_sample_id=rep(project.ID, nrow(voi_df)), voi_df, 
                         CL_ontology_id=gene.df$ontologyID[match(voi_df$cell_type_assignment,gene.df$cellName)])
  write.table(final.df, sep = "\t", quote = F, row.names = F,
              file = file.path(out_loc,"results/submission_table",paste0(ind.lib,"_metadata.tsv")))
  
  ## plotting the variables
  plot.list <- list()
  for (plot.type in voi){
    if (plot.type == "sctype_classification"){
      clrs <- ct.colors
    } else{
      clrs <- NULL
    }
    tryCatch({
      plot.list[[plot.type]] <- DimPlot(seu, reduction = "Xumap_", group.by = plot.type, 
                                        label = T, cols = clrs, repel = T) + 
        ggtitle(changeName.voi[match(plot.type, voi)])
    }, error=function(e){})
  }
  cowplot::plot_grid(plotlist = plot.list, nrow = n.row) + patchwork::plot_annotation(title = ind.lib) &  
    theme(plot.title = element_text(hjust = 0.5, size = 18, face="bold")) 
  ggsave(file.path(out_loc,"results/submission_table",paste0("multipanels_",ind.lib,".png")), width = 12, height = 5, bg = "white", dpi = 150)
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)

gene.df <- read.table(file.path(out_loc, "Azimuth_BM_level1.csv"), sep = ",", header = T)
ct_color <- c("darkorchid","skyblue2","dodgerblue2","gold","beige","sienna1","green4","navy",
              "chocolate4","red","darkred","#6A3D9A","maroon","yellow4","grey35","black","lightpink","grey80")
names(ct_color) <- c("B","CD4 T","CD8 T","DC","HSPC","Mono","NK","Other T","Macrophage",
                     "Early Eryth","Late Eryth","Plasma","Platelet","Stromal","Blast","Cancer","Pre Eryth","Unknown")

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id
purrr::walk(libraryID, ~ writeout(ind.lib = .x))
