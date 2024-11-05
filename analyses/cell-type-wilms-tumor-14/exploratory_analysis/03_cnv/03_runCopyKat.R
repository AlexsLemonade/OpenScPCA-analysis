library(dplyr)
library(Seurat)
library(copykat)
library(ggplot2)

path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 

library_id <- "SCPCL000850"

# create output dirs
scratch_out_dir <- file.path(path_anal, "scratch", "03_cnv","copykat")
dir.create(scratch_out_dir, showWarnings = FALSE, recursive = TRUE)
library_out_dir <- file.path(scratch_out_dir, library_id)
dir.create(library_out_dir, showWarnings = FALSE, recursive = TRUE)

# results_out_dir <- file.path(path_anal, "results", "03_cnv")
# dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)
# plots_out_dir <- file.path(path_anal, "plots", "03_cnv")
# dir.create(results_out_dir, showWarnings = FALSE, recursive = TRUE)


# load pre-processed sample objs & anchor transfer results
obj <- SeuratObject::LoadSeuratRds( file.path(path_anal,"scratch","00_preprocessing_rds",paste0(library_id,".rdsSeurat")) )
level <- "compartment"
predictions <- read.csv( file.path(path_anal, "results", "01_anchor_transfer_seurat", "RNA", paste0(library_id, "_", level,".csv")) ) 
obj <- AddMetaData(object = obj, metadata = predictions)

# prepare copykat input matrix & normal cell list
count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "count")
normal_cells <- predictions %>%
  tibble::column_to_rownames(var = "X") %>%
  filter(predicted.id == "immune")

# save copykat intermediate results by setting directory
setwd(library_out_dir)

# run copykat with reference normal cells
copykat_result <- copykat(
  rawmat = count_mat,
  id.type = "E",
  sam.name = library_id,
  norm.cell.names = rownames(normal_cells),
  ngene.chr = 2,
  plot.genes = FALSE,
  output.seg = FALSE,
  n.cores = 8
)
readr::write_rds(copykat_result, file = file.path(library_out_dir, paste0(library_id, "_copykat_resultobj.rds")))

# run copykat without normal cells
copykat_result_noref <- copykat(
  rawmat = count_mat,
  id.type = "E",
  sam.name = paste0(library_id,"_noref"),
  norm.cell.names = NULL,
  ngene.chr = 2,
  plot.genes = FALSE,
  output.seg = FALSE,
  n.cores = 8
)
readr::write_rds(copykat_result_noref, file = file.path(library_out_dir, paste0(library_id, "_noref_copykat_resultobj.rds")))

