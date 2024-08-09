# This script is to run the 01-clustering.Rmd template notebook over all samples

## Load libraries

library(SingleCellExperiment)
library(Seurat)      # to use Seurat workflow

library(tidyr)
library(dplyr)
library(DT)          # for table visualization

library(ggplot2)     # for visualization
library(patchwork)   # for visualization
library(ggplotify)
library(SCpubr)      # for visualization
library(assertthat)  # for visualization, to be able to use do_DimPlot
library(viridis)     # for visualization - colors

library(edgeR)       # For pseudobulk based differential expression (DE) analysis
library(DElegate)    # for pseudobulk based DE analysis and identification of marker genes

library(Azimuth)     # For label transfer

library(msigdbr)     # for gene set enrichment analysis or pathway enrichment
library(clusterProfiler)

library(data.table)

## define some config set-up / threshold

# We store in a config list names "cfg" parameters used all along the analysis to filter for p-value, log fold change, percentage of expression, etc.

cfg = list()
cfg$padj_thershold = 0.05
cfg$lfc_threshold = 1
cfg$rate1_threshold = 0.5
set.seed(12345)


## Define the base directories

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)

# The current data directory, found within the repository base directory
data_dir <- file.path(repository_base, "data", "2024-07-08", "SCPCP000006")

# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")


## Input and output files

#For this analysis, we worked with the _processed.rds data. 
#We builded a Seurat object based on the counts data and re-perform the analysis 
#[normalization --> reduction --> clustering] following the Seurat workflow.

# input
# use download-data.py to download the data

# get all files in the data directory
filelist <- list.files(data_dir, full.names = TRUE)

# select the 40 Wilms tumor single nucleus folders only
filelist <- filelist[! grepl(".tsv", filelist)]

# output
## html report for each of the Wilms tumor sample will be saved in 
path_to_report <- file.path(module_base, "notebook")
## some plots from the report will be saved in 
path_to_plot <- file.path(module_base, "/plots")
## after final decision on the clustering strategy, processed RDS files will be saved in 
path_to_result <- file.path(module_base, "/results")


# open the set of reference marker genes

CellType_metadata <- read.table("marker-sets/CellType_metadata.csv", header = TRUE, sep = ",")
Pathway_metadata <- read.csv("marker-sets/Pathway_metadata.csv", header = TRUE, sep = ",")


for (i in 2:length(filelist)) {
  
  
  rmarkdown::render(input = "notebook-template/01-clustering.Rmd", 
                    output_format = "html_document",
                    output_file = paste0("01-clustering_",gsub(paste0(data_dir, "/"), "", filelist[i]), ".html"),
                    output_dir = "notebook/01-clustering")
  
  
}

