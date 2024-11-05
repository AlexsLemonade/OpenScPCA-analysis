#!/usr/bin/env Rscript

#This script converts the "Azimuth_BM_level1.csv" into the submissio format for marker genes

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")

all_sources <- c("https://azimuth.hubmapconsortium.org/references/#Human%20-%20Bone%20Marrow",
                 "doi:10.1038/s41598-023-39152-z","https://sctype.app/database.php")
gene.df <- read.table(file.path(out_loc,"Azimuth_BM_level1.csv"), sep = ",", header = T)

df1 <- c()
for (i in 1:nrow(gene.df)){
  tmp.gene <- strsplit(gene.df$ensembl_id_positive_marker[i],",")[[1]]
  if (gene.df$cellName[i] == "Blast"){
    source <- all_sources[2]
  }else if (gene.df$cellName[i] %in% c("Pre Eryth","Cancer")){
    source <- all_sources[3]
  }else{source <- all_sources[1]}
  
  df1 <- rbind(df1,data.frame(ensembl_gene_id=tmp.gene,
                              cell_type=rep(gene.df$cellName[i],length(tmp.gene)),
                              source=rep(source,length(tmp.gene))))
}

df2 <- by(df1, df1$ensembl_gene_id, \(x) list2DF(lapply(x, \(.) toString(unique(.))))) |>
  do.call(what=rbind)

write.table(df2, file = file.path(out_loc,"submission_markerGenes.tsv"), sep = "\t", 
            row.names = F, quote = F)
