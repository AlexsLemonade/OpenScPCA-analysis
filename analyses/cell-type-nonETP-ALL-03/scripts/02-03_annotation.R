#!/usr/bin/env Rscript

#This script uses ScType to perform cell type annotation and CopyKat for tumor cells identification.
#We provides Azimuth_BM_level1.csv as the marker files for ScType, and generates sctype_classification metadata for the seurat object
#We also provides the top 10 celltypes with the highest sctype score for each cluster in the results/.
#Based on the annotated B cells, we provide them as the normal cells for CopyKat to identify tumor cells.
#The outputs are final rds file and two umap plots showing cell type and tumor/normal cells status.

library(dplyr)
library(Seurat)
library(HGNChelper)
library(copykat)
library(ggplot2)

# load gene set preparation function
#source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://github.com/IanevskiAleksandr/sc-type/raw/6db9eef49f185cf4d79bfec92a20fcf1edcccafb/R/sctype_score_.R")

gene_sets_prepare <- function(path_to_db_file, cell_type){
  cell_markers = read.csv(path_to_db_file, header = T)
  cell_markers = cell_markers[cell_markers$tissueType == cell_type,] 
  cell_markers$ensembl_id_positive_marker = gsub(" ","",cell_markers$ensembl_id_positive_marker); cell_markers$ensembl_id_negative_marker = gsub(" ","",cell_markers$ensembl_id_negative_marker)
  
  # correct gene symbols from the given DB (up-genes)
  cell_markers$ensembl_id_positive_marker = sapply(1:nrow(cell_markers), function(i){
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$ensembl_id_positive_marker[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(markers_all))}) #since the markers are provided in Ensembl ID, I removed checkGeneSymbols function here 
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  # correct gene symbols from the given DB (down-genes)
  cell_markers$ensembl_id_negative_marker = sapply(1:nrow(cell_markers), function(i){
    markers_all = gsub(" ", "", unlist(strsplit(cell_markers$ensembl_id_negative_marker[i],",")))
    markers_all = toupper(markers_all[markers_all != "NA" & markers_all != ""])
    markers_all = sort(markers_all)
    
    if(length(markers_all) > 0){
      suppressMessages({markers_all = unique(na.omit(markers_all))}) #since the markers are provided in Ensembl ID, I removed checkGeneSymbols function here 
      paste0(markers_all, collapse=",")
    } else {
      ""
    }
  })
  
  cell_markers$ensembl_id_positive_marker = gsub("///",",",cell_markers$ensembl_id_positive_marker);cell_markers$ensembl_id_positive_marker = gsub(" ","",cell_markers$ensembl_id_positive_marker)
  cell_markers$ensembl_id_negative_marker = gsub("///",",",cell_markers$ensembl_id_negative_marker);cell_markers$ensembl_id_negative_marker = gsub(" ","",cell_markers$ensembl_id_negative_marker)
  
  gs = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$ensembl_id_positive_marker[j]),",")))); names(gs) = cell_markers$cellName
  gs2 = lapply(1:nrow(cell_markers), function(j) gsub(" ","",unlist(strsplit(toString(cell_markers$ensembl_id_negative_marker[j]),",")))); names(gs2) = cell_markers$cellName
  
  list(gs_positive = gs, gs_negative = gs2)
}

running_scType <- function(gs_list, annot.obj, assay = "RNA", thres = 4){
  # check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
  seurat_package_v5 <- isFALSE('counts' %in% names(attributes(annot.obj[[assay]])));
  print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))
  # extract scaled scRNA-seq matrix
  scRNAseqData_scaled <- if (seurat_package_v5) as.matrix(annot.obj[[assay]]$scale.data) else as.matrix(annot.obj[[assay]]@scale.data)
  # run ScType (es.max is a table with #cell types as row and #cells as column, recording the sctype score for each cell/celltype combination)
  es.max <- sctype_score(scRNAseqData = scRNAseqData_scaled, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  # merge by cluster [repetitively finding top 10 cell types with the highest total ScType score for each cluster]
  ## es.max.cl - stores the descending value of total ScType score for each cell type in a particular cluster
  ## The following line prints out the top 10 cell types with the highest total ScType score for a cluster, along with #cells in that cluster
  ## which is then stored in cL_results (_sctype_top10_celltypes_perCluster.txt) - table with 4 columns (clusterID, cell type, total ScType score, #cells) 
  cL_results <- do.call("rbind", lapply(unique(annot.obj@meta.data$leiden_clusters), function(cl){
    es.max.cl = sort(rowSums(es.max[ ,rownames(annot.obj@meta.data[annot.obj@meta.data$leiden_clusters==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(annot.obj@meta.data$leiden_clusters==cl)), 10)
  }))
  sctype_scores <- cL_results %>% group_by(cluster) %>% top_n(n = 1, wt = scores)
  # set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/thres] <- "Unknown"
  print(sctype_scores[,1:3])
  
  annot.obj@meta.data$sctype_classification = ""
  for(cluster_num in unique(sctype_scores$cluster)){
    cl_type = sctype_scores[sctype_scores$cluster==cluster_num,];
    annot.obj@meta.data$sctype_classification[annot.obj@meta.data$leiden_clusters == cluster_num] = as.character(cl_type$type[1])
  }
  return (list(annot.obj, cL_results))
}

plot_modulescore <- function(gs_list, seu, sample.name){
  for (i in 1:length(gs_list$gs_positive)){
    seu <- AddModuleScore(object = seu, name = paste0(gsub(" ","",names(gs_list$gs_positive[i])),"_Features"), 
                          features = list(gs_list$gs_positive[[i]]))
  }
  features.vec <- paste0(gsub(" ","",names(gs_list$gs_positive)),"_Features1")
  Idents(seu) <- seu$leiden_clusters
  DotPlot(seu, features = features.vec) + theme(axis.text.x = element_text(angle=45,hjust=1))  +
    ggtitle(paste0(sample.name,": cell type module score")) +
    scale_color_gradient2(low = scales::muted("blue"), mid = "whitesmoke", high = scales::muted("red"), midpoint = 0)
  ggsave(file.path(out_loc,"plots",paste0(sample.name,"_features_dotplot.png")), width = 7, height = 7)
}

run_annot <- function(ind.lib){
  seu <- readRDS(file.path(out_loc,"results/rds",paste0(ind.lib,".rds")))
  final.obj <- seu #so I can keep an object without scale.data (it is taking too much space to save scale.data)
  seu <- ScaleData(seu, features = rownames(seu))
  gs_list <- gene_sets_prepare(db, tissue) #prepare gene sets
  ## relaxing the threshold of sctype_score larger than %ncells in a cluster, from 25% (default) to 10%
  res <- running_scType(gs_list, seu, thres = 10)   #cell type annotation
  seu$lowConfidence_annot <- res[[1]]$sctype_classification #seurat object with 'sctype_classification' metadata added
  res <- running_scType(gs_list, seu)
  seu <- res[[1]]
  #res[[2]] - A table with top 10 cell types with the highest scores for each cluster
  write.table(res[[2]], file = file.path(out_loc,"results",paste0(ind.lib,"_sctype_top10_celltypes_perCluster.txt")),
              row.names = F, sep = "\t", quote = F)
  p1 <- DimPlot(seu, reduction = "Xumap_", group.by = "sctype_classification", cols = ct_color) +
    ggtitle(paste0(ind.lib,": cell type"))
  p2 <- DimPlot(seu, reduction = "Xumap_", group.by = "lowConfidence_annot", cols = ct_color) +
    ggtitle(paste0(ind.lib,": low confidence annotation"))
  p1 + p2
  ggsave(file.path(out_loc,"plots",paste0(ind.lib,"_celltype.png")), width = 10, height = 5)

  plot_modulescore(gs_list, seu, ind.lib)
  
  #using copykat for tumor cells identification
  norm.cells <- colnames(seu)[which(seu$sctype_classification=="B")]
  n_cores <- parallel::detectCores() - 1
  if (length(norm.cells) > 0){    #the sample has B cells annotated
    copykat.test <- copykat(rawmat=seu@assays[["RNA"]]@counts, id.type="Ensemble",
                            ngene.chr=5, win.size=25, KS.cut=0.1, sam.name=ind.lib,
                            distance="euclidean", norm.cell.names=norm.cells,
                            output.seg="FALSE", plot.genes="TRUE", genome="hg20",n.cores=n_cores)
    idx <- match(colnames(seu), copykat.test$prediction$cell.names)
    seu$copykat.pred <- copykat.test$prediction$copykat.pred[idx]
    DimPlot(seu, reduction = "Xumap_", group.by = "copykat.pred") +
      ggtitle(paste0(ind.lib,": copykat prediction"))
    ggsave(file.path(out_loc,"plots",paste0(ind.lib,"_copykatPred.png")), width = 7, height = 7)
    
    voi <- c('leiden_clusters','sctype_classification','lowConfidence_annot','copykat.pred')
    final.obj$copykat.pred <- seu$copykat.pred
  }else{
    voi <- c('leiden_clusters','sctype_classification','lowConfidence_annot')
  }
  
  write.table(data.frame(FetchData(seu, vars = voi)), sep = "\t", quote = F,
              file = file.path(out_loc,"results",paste0(ind.lib,"_metadata.txt")))
  final.obj$sctype_classification <- seu$sctype_classification
  final.obj$lowConfidence_annot <- seu$lowConfidence_annot
  saveRDS(final.obj, file = file.path(out_loc,"results/rds",paste0(ind.lib,".rds")))
}


project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id
# DB file
db <- file.path(out_loc,"Azimuth_BM_level1.csv")
tissue <- "Immune system"  
ct_color <- c("darkorchid","skyblue2","dodgerblue2","gold","beige","sienna1","green4","navy",
              "chocolate4","red","darkred","#6A3D9A","maroon","yellow4","grey35","black","lightpink","grey80")
names(ct_color) <- c("B","CD4 T","CD8 T","DC","HSPC","Mono","NK","Other T","Macrophage",
                     "Early Eryth","Late Eryth","Plasma","Platelet","Stromal","Blast","Cancer","Pre Eryth","Unknown")

purrr::walk(libraryID, run_annot)
