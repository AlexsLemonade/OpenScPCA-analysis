#!/usr/bin/env Rscript

# USAGE:
# Rscript 00_run_workflow.R
#
# USAGE in CI:
# Rscript 00_run_workflow.R --testing

# Set up options ----------------------------------------------------------------

library(optparse)

option_list <- list(
  make_option(
    opt_str = c("--testing"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Use this flag when running on test data"
  )
)

opts <- parse_args(OptionParser(option_list = option_list))
running_ci <- opts$testing

# Run the Label transfer from two fetal references ------------------------------

# get list of samples in the library --------------------------------------------
root_dir <- rprojroot::find_root(rprojroot::is_git_root)
project_id <- "SCPCP000006"
sample_metadata_file <- file.path(root_dir, "data", "current", project_id, "single_cell_metadata.tsv")
metadata <- read.table(sample_metadata_file, sep = "\t", header = TRUE)

# set path to this module--------------------------------------------------------
module_base <- file.path(root_dir, "analyses", "cell-type-wilms-tumor-06")

# Download and create the fetal kidney reference (Stewart et al) ----------
system(command = glue::glue("Rscript ", file.path(module_base,"scripts", "download-and-create-fetal-kidney-ref.R")))


# Characterize the two fetal references -----------------------------------------

# Characterize the fetal full reference (Cao et al.)
# To be done, next PR
notebook_template_dir <- file.path(module_base, "notebook_template")
notebook_output_dir <- file.path(module_base, "notebook")
# Characterize the fetal kidney reference (Stewart et al.)
rmarkdown::render(input = file.path(notebook_template_dir, "00b_characterize_fetal_kidney_reference_Stewart.Rmd"),
                  output_format = "html_document",
                  output_file = "00b_characterization_fetal_kidney_reference_Stewart.html",
                  output_dir = file.path(notebook_output_dir, "00-reference"))


# Run the workflow for (all) samples in the project -----------------------------
for (sample_id in metadata$scpca_sample_id) {

  # create a directory to save the pre-processed and labeled `Seurat` objects
  dir.create(file.path(module_base, "results", sample_id), showWarnings = FALSE)
  # create a directory to save the notebooks
  dir.create(file.path(module_base, "notebook", sample_id), showWarnings = FALSE)

  # Pre-process the data - `Seurat` workflow
  rmarkdown::render(input = file.path(notebook_template_dir, "01_seurat-processing.Rmd"),
                    params = list(scpca_project_id = project_id, sample_id = sample_id),
                    output_format = "html_document",
                    output_file = paste0("01_seurat_processing_", sample_id, ".html"),
                    output_dir = file.path(notebook_output_dir,  sample_id))

  if (!running_ci) {
    # Label transfer from the Cao reference using Azimuth
    rmarkdown::render(input = file.path(notebook_template_dir, "02a_label-transfer_fetal_full_reference_Cao.Rmd"),
                      params = list(scpca_project_id = project_id, sample_id = sample_id),
                      output_format = "html_document",
                      output_file = paste0("02a_fetal_all_reference_Cao_", sample_id, ".html"),
                      output_dir = file.path(notebook_output_dir, sample_id))

    # Label transfer from the Stewart reference using Seurat
    rmarkdown::render(input = file.path(notebook_template_dir, "02b_label-transfer_fetal_kidney_reference_Stewart.Rmd"),
                      params = list(scpca_project_id = project_id, sample_id = sample_id),
                      output_format = "html_document",
                      output_file = paste0("02b_fetal_kidney_reference_Stewart_", sample_id, ".html"),
                      output_dir = file.path(notebook_output_dir, sample_id))
    
    # Cluster exploration
    rmarkdown::render(input = file.path(notebook_template_dir, "03_clustering_exploration.Rmd"),
                      params = list(scpca_project_id = project_id, sample_id = sample_id),
                      output_format = "html_document",
                      output_file = paste0("03_clustering_exploration_", sample_id, ".html"),
                      output_dir = file.path(notebook_output_dir, sample_id))

  }
}

if (!running_ci) {
  rmarkdown::render(input = file.path(notebook_output_dir, "04_annotation_Across_Samples_exploration.Rmd"),
                    output_format = "html_document",
                    output_file = "04_annotation_Across_Samples_exploration.html",
                    output_dir = notebook_output_dir)
}


# We run copyKAT for a subselection of samples selected in "04_annotation_Across_Samples_exploration.Rmd" with and without healthy cells as reference

for (sample_id in c("SCPCS000179",
                    "SCPCS000184",
                    "SCPCS000194", 
                    "SCPCS000205",
                    "SCPCS000208")){
  
  # We run and explore copykat using euclidian distance parameter
  system(command = glue::glue("Rscript ", file.path(module_base,"scripts", "05_copyKAT.R"), " --sample_id ", sample_id, " --n_core ", 32, " --distance ", "euclidean"))
  
  rmarkdown::render(input = file.path(notebook_template_dir, "05_cnv_copykat_exploration.Rmd"),
                    params = list(sample_id = sample_id, distance = "euclidean"),
                    output_format = "html_document",
                    output_file = paste0("05_cnv_copykat_euclidean_exploration_", sample_id, ".html"),
                    output_dir = file.path(notebook_output_dir, sample_id))
  
  # We run and explore copykat using spearman distance parameter
  
  system(command = glue::glue("Rscript ", file.path(module_base,"scripts", "05_copyKAT.R"), " --sample_id ", sample_id, " --n_core ", 32, " --distance ", "spearman"))
  
  rmarkdown::render(input = file.path(notebook_template_dir, "05_cnv_copykat_exploration.Rmd"),
                    params = list(sample_id = sample_id, distance = "spearman"),
                    output_format = "html_document",
                    output_file = paste0("05_cnv_copykat_spearman_exploration_", sample_id, ".html"),
                    output_dir = file.path(notebook_output_dir, sample_id))
  
  # run infercnv with different selection of normal cells as reference
  for (sample_id in c(
                      "SCPCS000184",
                      "SCPCS000194", 
                      "SCPCS000205",
                      "SCPCS000208")){

  # We run and explore infercnv using immune cells as reference
  system(command = glue::glue("Rscript ", file.path(module_base,"scripts", "06_infercnv.R"), " --sample_id ", sample_id, " --reference ", "immune"))
  

   # We run and explore infercnv using endothelial cells as reference
  system(command = glue::glue("Rscript ", file.path(module_base,"scripts", "06_infercnv.R"), " --sample_id ", sample_id, " --reference ", "endothelium"))
  
  # We run and explore infercnv using immune and endothelium cells as reference
  system(command = glue::glue("Rscript ", file.path(module_base,"scripts", "06_infercnv.R"), " --sample_id ", sample_id, " --reference ", "both"))
  
  # We run and explore infercnv using no normal reference
  system(command = glue::glue("Rscript ", file.path(module_base,"scripts", "06_infercnv.R"), " --sample_id ", sample_id, " --reference ", "none"))
  
  
  rmarkdown::render(input = file.path(notebook_template_dir, "06_cnv_infercnv_exploration.Rmd"),
                    params = list(sample_id = sample_id),
                    output_format = "html_document",
                    output_file = paste0("06_cnv_exploration_", sample_id, ".html"),
                    output_dir = file.path(notebook_output_dir, sample_id))
  
  
}

