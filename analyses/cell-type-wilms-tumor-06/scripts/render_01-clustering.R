#!/usr/bin/env Rscript

# get list of samples in the library
root_dir <- rprojroot::find_root(rprojroot::is_git_root)
sample_metadata_file <- file.path(root_dir, "data", "current", "SCPCP000006", "single_cell_metadata.tsv")
metadata <- read.table(sample_metadata_file, sep = "\t", header = TRUE)

for (i in metadata$scpca_sample_id[9:11]) {
  
  rmarkdown::render(input = "notebook-template/01-clustering.Rmd", 
                    params = list(scpca_project_id = metadata$scpca_project_id[metadata$scpca_sample_id ==i], sample_id = i),
                    output_format = "html_document",
                    output_file = paste0("01-clustering_",i, ".html"),
                    output_dir = "notebook/01-clustering")
  
  
}

