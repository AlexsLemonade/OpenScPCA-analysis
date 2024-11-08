#!/usr/bin/env Rscript

## Calculates silhouette score and purity for each cluster and evaluates their stability,
## using the functions available in evaluate-cluster.R (on OpenscPCA portal)

library(Seurat)
library(dplyr)

run_eval <- function(ind.lib) {
  seu <- readRDS(file.path(out_loc, "results/rds", paste0(ind.lib, ".rds")))
  clusID.df <- data.frame(FetchData(seu, vars = "leiden_clusters")) |> tibble::rownames_to_column(var = "cell_id")
  colnames(clusID.df)[2] <- "cluster"
  cluster_df1 <- rOpenScPCA::calculate_silhouette(x = seu, cluster_df = clusID.df, pc_name = "Xpca_")
  cluster_df2 <- rOpenScPCA::calculate_purity(x = seu, cluster_df = clusID.df, pc_name = "Xpca_")
  final_df <- merge(cluster_df1, cluster_df2, by = c("cell_id", "cluster"))
  perClus_df <- final_df %>%
    group_by(cluster) %>%
    summarise(avgSil = mean(silhouette_width), avgPur = mean(purity)) %>%
    data.frame()
  stability_df <- rOpenScPCA::calculate_stability(
    x = seu, clusters = clusID.df,
    pc_name = "Xpca_", algorithm = "leiden",
    resolution = 1.0, objective_function = "modularity",
    seed = 10
  )
  write.table(final_df,
    sep = "\t", row.names = F, quote = F,
    file = file.path(out_loc, "results/evalClus/", paste0(ind.lib, "_sil-purity_perCell.txt"))
  )
  write.table(stability_df,
    sep = "\t", row.names = F, quote = F,
    file = file.path(out_loc, "results/evalClus/", paste0(ind.lib, "_stability.txt"))
  )
  write.table(perClus_df,
    sep = "\t", row.names = F, quote = F,
    file = file.path(out_loc, "results/evalClus/", paste0(ind.lib, "_avgSil-purity_perClus.txt"))
  )
}

project_root <- rprojroot::find_root(rprojroot::is_git_root)
projectID <- "SCPCP000003"
out_loc <- file.path(project_root, "analyses/cell-type-ETP-ALL-03")
data_loc <- file.path(project_root, "data/current", projectID)
dir.create(file.path(out_loc, "results/evalClus"), showWarnings = FALSE)

metadata <- read.table(file.path(data_loc, "single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
  metadata$diagnosis == "Early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id

purrr::walk(libraryID, run_eval)
