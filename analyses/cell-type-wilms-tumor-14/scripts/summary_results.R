#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(purrr)
library(tibble)

option_list <- list(
  make_option(
    opt_str = c("--metadata"),
    type = "character",
    default = NULL,
    help = "Path to cohort metadata"
  ),
  make_option(
    opt_str = c("--testing"),
    type = "logical",
    default = FALSE,
    action = "store_true",
    help = "Use this flag when running on test data"
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))
running_ci <- opt$testing

if (running_ci) {
  assays <- c("RNA")
} else {
  assays <- c("RNA", "SCT")
}

path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 
# path_meta <- file.path(path_repo,"data","current","SCPCP000014","single_cell_metadata.tsv") # keep for debug
path_meta <- file.path(opt$metadata)
meta <- read.table(path_meta, sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

# output table
outfile <- file.path(file.path(path_anal, "results", "01_anchor_transfer_seurat", "summary_results.csv"))

# for testing
# anchor_assay <- "RNA"; level <- "compartment"; 
# library_id <- "SCPCL000846"; sample_id <- "SCPCS000514"

################### Functions ################### 
read_library_csv <- function(library_id, 
                             sample_id,
                             anchor_assay,
                             level) {
  # get results from anchor transfer
  anchor_dir <- file.path(path_anal, "results", "01_anchor_transfer_seurat", anchor_assay)
  col_name <- paste("cell_type_assignment",anchor_assay,level,sep = "_")
  df <- read.csv( file.path(anchor_dir, paste0(library_id, "_", level, ".csv")) ) %>%
    select(c("X","predicted.id"))
  colnames(df) <- c("cell_barcode", col_name)
  df <- df %>%
    tibble::add_column(scpca_sample_id = sample_id, scpca_library_id = library_id, .before = col_name )
  return(df)
}


concat_sample <- function(anchor_assay,
                         level,
                         meta) {
  prediction_concat <- purrr::map2(
    meta$scpca_library_id,
    meta$scpca_sample_id,
    \(library_id, sample_id) read_library_csv(library_id = library_id,
                                              sample_id = sample_id,
                                              anchor_assay = anchor_assay,
                                              level = level)
  )
  prediction_concat <- prediction_concat %>%
    bind_rows()
  return(prediction_concat)
}


################### Creating summary table ################### 


arg_df <- tidyr::expand_grid(
  anchor_assay = assays,
  level = c("compartment", "celltype")
)
dfs <- purrr::pmap(arg_df,
             \(anchor_assay, level) concat_sample(anchor_assay, level, meta)
)
# add rough tumor cell classification, CM = tumor cells, immune/endothelium = normal
normal_compartments <- c("immune", "endothelium")
final_table <- dfs %>%
  purrr::reduce(left_join) %>%
  mutate( tumor_cell_classification_RNA = case_when(
    cell_type_assignment_RNA_compartment %in% normal_compartments ~ "normal",
    cell_type_assignment_RNA_compartment == "fetal_nephron" & cell_type_assignment_RNA_celltype == "Cap mesenchyme" ~ "tumor",
    .default = "unknown"
  ) )
if ("SCT" %in% assays) {
  final_table <- final_table %>%
    mutate( tumor_cell_classification_SCT = case_when(
      cell_type_assignment_SCT_compartment %in% normal_compartments ~ "normal",
      cell_type_assignment_SCT_compartment == "fetal_nephron" & cell_type_assignment_SCT_celltype == "Cap mesenchyme" ~ "tumor",
      .default = "unknown"
    ) )
}
final_table <- final_table %>%
  # use logNormalize as "default" final result, combine immune, endo, and stroma compartments
  mutate( cell_type_assignment = if_else(cell_type_assignment_RNA_compartment != "fetal_nephron", 
                                         cell_type_assignment_RNA_compartment, 
                                         cell_type_assignment_RNA_celltype ) ) %>%
  mutate( tumor_cell_classification = tumor_cell_classification_RNA )


write.csv(final_table, file = outfile)
