library(dplyr)
library(Seurat)
library(infercnv)
library(ggplot2)

path_repo <- rprojroot::find_root(rprojroot::is_git_root)
path_anal <- file.path(path_repo,"analyses","cell-type-wilms-tumor-14") 

library_id <- "SCPCL000850"

# create output dirs
scratch_out_dir <- file.path(path_anal, "scratch", "03_cnv","infercnv")
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

# create annotation file for infercnv
annotation_file <- file.path(scratch_out_dir, paste0(library_id,"_annotations_file_infercnv.txt"))
# annotation_df <- data.frame(cell = obj$barcodes, cluster = obj$cluster)
annotation_df <- data.frame(cell = obj$barcodes, cluster = obj$predicted.id)
write.table(annotation_df, file = annotation_file, sep = "\t", quote = F, row.names = F, col.names = F)
# check gene order file, from infercnv ftp
gene_order_file <- file.path(scratch_out_dir,"Homo_sapiens.GRCh38.104.gene_order.txt")
gene_order <- read.table(file = gene_order_file) %>%
  arrange(V2, V3)

obj <- obj[rownames(obj) %in% gene_order$V1]
count_mat <- SeuratObject::GetAssayData(obj, assay = "RNA", layer = "count")

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=count_mat,
                                    annotations_file=annotation_file,
                                    delim="\t",
                                    gene_order_file=gene_order_file,
                                    ref_group_names="immune") 


infercnv_out = infercnv::run(infercnv_obj,
                             cutoff = 0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir = library_out_dir, 
                             denoise = TRUE,
                             HMM = TRUE,
                             save_rds = FALSE,
                             num_threads = 8)

readr::write_rds(infercnv_out, file = file.path(library_out_dir,"infercnv_out.rds"))