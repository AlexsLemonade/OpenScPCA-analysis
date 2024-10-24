#!/usr/bin/env Rscript

#sudo apt-get install r-cran-rjags
#This script runs inferCNV using the new B cells, and trying to identify "non-malignant" group based on the infercnv.png

library(Seurat)

run_inferCNV <- function(ind.lib){
  dir.create(file.path(scratch_dir, ind.lib), showWarnings = FALSE)
  seu <- readRDS(file.path(out_loc,"results/rds",paste0(ind.lib,".rds")))
  annot.file <- file.path(out_loc,"results",paste0(ind.lib,"_newB-normal-annotation.txt"))
  infercnv_obj <- infercnv::CreateInfercnvObject(raw_counts_matrix=seu@assays[["RNA"]]@counts,
                                                 annotations_file=annot.file,
                                                 delim="\t",
                                                 gene_order_file=geneFile,
                                                 ref_group_names="new B") 
  options(scipen = 100)
  infercnv_obj <- infercnv::run(infercnv_obj,
                                cutoff = 0.1, # use 1 for smart-seq, 0.1 for 10x-genomics
                                out_dir = file.path(scratch_dir,ind.lib), # save all intermediate files to scratch dir
                                denoise = T, HMM = T, analysis_mode = "samples",
                                save_rds = F, # don't save the intermediate rds files
                                num_threads = parallel::detectCores() - 1)
  
  # create table with barcodes and CNVs for each chromosome
  infercnv::add_to_seurat(
    seurat_obj = NULL,
    infercnv_output_path = file.path(scratch_dir, ind.lib))
  
  ### adding clusterID from the cutting of hierarchical clustering to seu object
  final_cnv_obj <- readRDS(file.path(out_loc, "results/infercnv_output", paste0(ind.lib,"_run.final.infercnv_obj")))
  hres <- final_cnv_obj@tumor_subclusters$hc$all_observations
  obs.clusID <- cutree(hres,4) #4 clusters for SCPCL000703
  hres <- dendextend::color_branches(hres, k= 4) 
  plot(hres)
  seu$infercnv.pred <- obs.clusID[match(names(obs.clusID), colnames(seu))]
}

project_root  <- rprojroot::find_root(rprojroot::is_git_root)
out_loc <- file.path(project_root, "analyses/cell-type-nonETP-ALL-03")
data_loc <- file.path(project_root, "data/current",projectID)
scratch_dir <- file.path(out_loc,"scratch")

metadata <- read.table(file.path(data_loc,"single_cell_metadata.tsv"), sep = "\t", header = T)
metadata <- metadata[which(metadata$scpca_project_id == projectID &
                             metadata$diagnosis == "Non-early T-cell precursor T-cell acute lymphoblastic leukemia"), ]
libraryID <- metadata$scpca_library_id

geneFile <- file.path(project_root, "references/infercnv_refs/Homo_sapiens.GRCh38.104.gene_order.txt")

purrr::walk(libraryID, run_inferCNV)
