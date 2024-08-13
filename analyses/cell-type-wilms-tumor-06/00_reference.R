# get list of samples in the library
metadata <- read.table("sample_metadata.tsv", sep = "\t", header = TRUE)


# Render the reports for (all) samples in the project
for (i in metadata$scpca_sample_id[9:11]) {
  # Label transfer from the Cao reference using Azimuth
  rmarkdown::render(input = "notebook_template/01_fetal_reference_Cao_Azimuth.Rmd", 
                    params = list(scpca_project_id = metadata$scpca_project_id[metadata$scpca_sample_id ==i], sample_id = i),
                    output_format = "html_document",
                    output_file = paste0("01_fetal_reference_Cao_Azimuth_",i, ".html"),
                    output_dir = "notebook/00-reference")
  
  # Label transfer from the Stewart reference using Seurat
  rmarkdown::render(input = "notebook_template/01_fetal_reference_Stewart_Seuratv4.Rmd", 
                    params = list(scpca_project_id = metadata$scpca_project_id[metadata$scpca_sample_id ==i],sample_id = i),
                    output_format = "html_document",
                    output_file = paste0("01_fetal_reference_Stewart_Seuratv4_",i, ".html"),
                    output_dir = "notebook/00-reference")
  
  # Label transfer from the Stewart reference using Azimuth
  # Note this report need the download of the azimuth compatible reference in results s3 bucket.
  # I would be in favor to remove this part in the module analysis as the result with Seurat v4 are similar.
  # I just left it here now for information
  rmarkdown::render(input = "notebook_template/01_fetal_reference_Stewart_Azimuthv5.Rmd", 
                    params = list(scpca_project_id = metadata$scpca_project_id[metadata$scpca_sample_id ==i], sample_id= i),
                    output_format = "html_document",
                    output_file = paste0("01_fetal_reference_Stewart_Azimuthv5_",i, ".html"),
                    output_dir = "notebook/00-reference")
  
  }



