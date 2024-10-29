# Define path -------------------------------------------------------------------
# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")

# Define threshold score. Only cells with annotations greater than this threshold will be used
score_threshold <- 0.85

# Select samples that haven't be pre-treated with chemotherapies
# Even if normal cells shouldn't be affected by chemotherapy in terns of CNV, we decided to only take endothelial and immune cells from non-treated samples to build the reference of normal cells

sample_upfront_resection <- readr::read_tsv(
  file.path( "..", "..", "data", "current", "SCPCP000006", "single_cell_metadata.tsv"),
  show_col_types = FALSE
) |>
  dplyr::filter(treatment == "Upfront resection") |>
  dplyr::pull(scpca_sample_id)

# For each upfront resection samples, which subset the `Seurat` object to immune and endothelial cells
s <- list()

for(sample_id in sample_upfront_resection){
  # define result directory
  result_dir               <- file.path(module_base, "results", sample_id)
  # load `Seurat` object
  srat_tmp                 <- readRDS(file.path(result_dir,  paste0("02b-fetal_kidney_label-transfer_",  sample_id, ".Rds")))
  # proceed only if we have more than one normal cell in the sample, else the merge won't work
  
    # turn the `SCT` slot to NULL to avoid ERROR when merging (https://github.com/satijalab/seurat/issues/6462#issuecomment-1384499533)
    srat_tmp[["SCT"]]        <- NULL
    # subset only immune and endothelium cells
    srat_tmp               <- subset(srat_tmp, subset = fetal_kidney_predicted.compartment %in% c("endothelium", "immune") )
      if(dim(srat_tmp)[2]>1){
      s[[sample_id]]         <- srat_tmp
      # rename cells to indicate that these cells are used as spike in of normal reference
      colnames(s[[sample_id]]) <- paste("spike", sample_id, colnames(s[[sample_id]]), sep = "_")
  }
  
}

# merge all normal cells in one `Seurat` object
srat_normal <- base::merge(s[[1]], y = s[2:length(s)])

# we only keep as nornal cells immune and endothelial cells with a high mapping score
srat_normal <- subset(srat_normal, subset = fetal_kidney_predicted.compartment.score > 0.85)
# save the `Seurat` object

srat_normal_file <- file.path(module_base, "results","references", '06b_normal-cell-reference.rds')
saveRDS(srat_normal, srat_normal_file)
