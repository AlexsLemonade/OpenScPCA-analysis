#!/usr/bin/env Rscript

# Run `infercnv` for one sample with or without a healthy reference
# infercnv
#
# USAGE:
# Rscript 06_infercnv.R \
#   --sample_id SCPCS000194
#   --ncore 16
#

library(optparse)
library(Seurat)
library(infercnv)
library(reshape2)
# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("-s", "--sample_id"),
    type = "character",
    default = "SCPCS000205",
    help = "The sample_id of the sample to be used for inference of genomic copy number using infercnv "
  ),
  make_option(
    opt_str = c("-r", "--reference"),
    type = "character",
    default = NULL,
    help = "Reference cells to use as normal cells, either none, immune, endothelium or both"
  ),
  make_option(
    opt_str = c("-g", "--gene_order_file"),
    type = "character",
    default = NULL,
    help = "Path to gene order file as tab delimited .txt file with no column headers.
      Columns are: Ensembl gene id, chr, start, stop."
  )
)

opts <- parse_args(OptionParser(option_list = option_list))

# paths to data ----------------------------------------------------------------

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")
# Path to the result directory
result_dir <- file.path(module_base, "results", opts$sample_id)
# path to output infercnv object

output_dir <- file.path(result_dir,  "06_infercnv", "reference-"opts$reference )



# Read in data -----------------------------------------------------------------
srat <- readRDS(
  file.path(result_dir,  paste0("02b-fetal_kidney_label-transfer_",  opts$sample_id, ".Rds"))
)

# Extract raw counts -----------------------------------------------------------
counts <- GetAssayData(object = srat, assay = "RNA", layer = "counts")

# Create a dataframe of annotation ---------------------------------------------
annot_df <- data.frame(condition = as.character(srat$fetal_kidney_predicted.compartment))
rownames(annot_df) <- colnames(counts)
stopifnot("Incorrect reference provided" = reference %in% c("none", "immune", "endothelium", "both")
if(opts$reference == "both"){
  normal_cells <- c("endothelium", "immune")
} else if(opts$reference == "none"){
  normal_cells <- NULL
} else{
  normal_cells <- opts$reference
}

# create output directory if it does not exist
dir.create(output_dir, recursive = TRUE)

# make sure we have a genome position file
if (is.null(opts$gene_order_file)){
  gene_order_file <- file.path("scratch", 'gencode_v19_gene_pos.txt')
} else{
  gene_order_file <- opts$gene_order_file
}
if (!file.exists(gene_order_file)) {
  download.file(url = 'https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v21_gen_pos.complete.txt',
                destfile = gene_order_file)
  gene_order <- read_tsv(gene_order_file, col_names =  FALSE)
  tmp <- gsub("\\..*","",gene_order$X1) 
  tmp <- gsub(".*\\|","",tmp) 
  gene_order$X1 <- tmp
  gene_order <- gene_order[grepl("ENSG", x = gene_order$X1),]
  
  write_tsv(col_names = FALSE, gene_order, gene_order_file, append = FALSE)
}

# Run infercnv ------------------------------------------------------------------
# create inferCNV object and run method

infercnv_obj = infercnv::CreateInfercnvObject(
  raw_counts_matrix = as.matrix(counts),
  annotations_file = annot_df,
  ref_group_names = normal_cells,
  gene_order_file = gene_order_file)

infercnv_obj = infercnv::run(
  infercnv_obj,
  cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
  out_dir=output_dir, 
  cluster_by_groups=T, 
  plot_steps=F,
  denoise=T,
  sd_amplifier=3,  # sets midpoint for logistic
  noise_logistic=TRUE,
  HMM = F  # turns gradient filtering on
)
