# This script downloads and adapt the genome position file for use with infercnv

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")

# load library
library(tidyverse)
library(dplyr)

# Define the gene order file
gene_order_file <- file.path(module_base, "results","references", 'gencode_v19_gene_pos.txt')

# Define the arm order file
arm_order_file <- file.path(module_base, "results","references", 'hg38_cytoBand.txt')

# Define the gene/arm combined information file
gene_arm_order_file <- file.path(module_base, "results","references", 'gencode_v29_gene_pos_arm.txt')

# Add to the gene order file information on the chromosome arms
  if (!file.exists(gene_arm_order_file)) {
    
    # Download the the gene order file
    if (!file.exists(gene_order_file)) {
      
      download.file(url = 'https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt',
                  destfile = gene_order_file)
    }
    gene_order <- read_tsv(gene_order_file, col_names =  FALSE)
    tmp <- gsub("\\..*","",gene_order$X1) 
    tmp <- gsub(".*\\|","",tmp) 
    gene_order$X1 <- tmp
    gene_order <- gene_order[grepl("ENSG", x = gene_order$X1),]
    
    # Download the the arm order file
    if (!file.exists(arm_order_file)) {
    download.file(url = 'http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz',
                  destfile = arm_order_file)
    }
    # Load cytoBand file into R
    cytoBand <- read.table(arm_order_file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
    
    # Assign column names
    colnames(cytoBand) <- c("chrom", "start", "end", "band", "stain")
    
    # Add a column for the chromosome arm (p or q)
    cytoBand <- cytoBand %>%
      mutate(arm = substr(band, 1, 1)) # Extract 'p' or 'q' from the band column
    
    # Group by chromosome and arm to calculate arm start and end positions
    chromosome_arms <- cytoBand %>%
      group_by(chrom, arm) %>%
      summarize(
        start = min(start),
        end = max(end),
        .groups = "drop"
      )
    
    # Filter out unwanted chromosomes (e.g., scaffolds, mitochondrial DNA)
    chromosome_arms <- chromosome_arms %>%
      filter(grepl("^chr[0-9XY]+$", chrom)) # Keep only standard chromosomes
    
   gene_order$X5 <- "p"
    for(i in 1:dim(gene_order)[1]){
      start_q <- chromosome_arms[chromosome_arms$chrom %in% gene_order[i, 2], ]
      if(dim(start_q)[1]>0){
       start_q <- start_q[start_q$arm == "p", 4]
        if(gene_order[i,3] > start_q){
        gene_order[i,5] <- "q"
        }
      }else{
        gene_order[i,5] <- "NA"
      }
      gene_order[i,2] <- paste0(gene_order[i,2], gene_order[i,5])
      
    }
   gene_order$X2 <- factor(gene_order$X2, levels = c("chr1p",
                                                     "chr1q", 
                                                     "chr2p",
                                                     "chr2q",
                                                     "chr3p",
                                                     "chr3q",
                                                     "chr4p",
                                                     "chr4q",
                                                     "chr5p",
                                                     "chr5q",
                                                     "chr6p",
                                                     "chr6q",
                                                     "chr7p",
                                                     "chr7q",
                                                     "chr8p",
                                                     "chr8q",
                                                     "chr9p",
                                                     "chr9q", 
                                                     "chr10p",
                                                     "chr10q", 
                                                     "chr11p",
                                                     "chr11q", 
                                                     "chr12p",
                                                     "chr12q",
                                                     "chr13p",
                                                     "chr13q",
                                                     "chr14p",
                                                     "chr14q",
                                                     "chr15p",
                                                     "chr15q",
                                                     "chr16p",
                                                     "chr16q",
                                                     "chr17p",
                                                     "chr17q",
                                                     "chr18p",
                                                     "chr18q",
                                                     "chr19p",
                                                     "chr19q",
                                                     "chr20p",
                                                     "chr20q",
                                                     "chr21p",
                                                     "chr21q",
                                                     "chr22p",
                                                     "chr22q",
                                                     "chrXp",
                                                     "chrXq",
                                                     "chrYp",
                                                     "chrYq",
                                                     "chrMNA"))
   gene_order <- gene_order[order(gene_order$X2, gene_order$X3),]
   
   write_tsv(col_names = FALSE, gene_order[,-5], gene_arm_order_file, append = FALSE)
   
  }
  