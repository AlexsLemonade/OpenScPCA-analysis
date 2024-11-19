#!/usr/bin/env Rscript

#This script prepares the gene set for each cell type, extracting them from the marker list 

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