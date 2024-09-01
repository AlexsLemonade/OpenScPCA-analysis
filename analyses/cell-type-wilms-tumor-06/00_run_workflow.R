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
  

    # Label transfer from the Stewart reference using Seurat
    rmarkdown::render(input = file.path(module_base, "notebook_template", "02b_label-transfer_fetal_kidney_reference_Stewart.Rmd"),
                      params = list(scpca_project_id = project_id, sample_id = sample_id),
                      output_format = "html_document",
                      output_file = paste0("02b_fetal_kidney_reference_Stewart_", sample_id, ".html"),
                      output_dir = file.path(notebook_output_dir, sample_id))

  }
}



