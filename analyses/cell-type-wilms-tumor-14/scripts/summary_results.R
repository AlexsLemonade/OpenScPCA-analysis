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
  )
)

# Parse options
opt <- parse_args(OptionParser(option_list = option_list))

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
                         level) {
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
  anchor_assay = c("RNA", "SCT"),
  level = c("compartment", "celltype")
)
dfs <- purrr::pmap(arg_df,
             \(anchor_assay, level) concat_sample(anchor_assay = anchor_assay, level = level)
)
final_table <- dfs %>%
  purrr::reduce(left_join)

write.csv(final_table, file = outfile)
