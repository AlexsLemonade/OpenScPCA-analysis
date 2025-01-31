# This script downloads and adapt the genome position file for use with infercnv

# load library
library(dplyr)
library(stringr)


# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")

# Define the file paths
gene_order_file <- file.path(module_base, "results", "references", "gencode_v38_gene_pos.txt")
arm_order_file <- file.path(module_base, "results", "references", "hg38_cytoBand.txt")
gene_arm_order_file <- file.path(module_base, "results", "references", "gencode_v38_gene_pos_arm.txt")

if (!file.exists(gene_arm_order_file)) {
  
  # Download the updated gene order file for GENCODE v38 (Ensembl 104)
  if (!file.exists(gene_order_file)) {
    download.file(url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.annotation.gtf.gz",
                  destfile = gene_order_file)
  }
  
  # Read the GTF file and extract gene positions
  gtf_data <- read.table(gzfile(gene_order_file), header = FALSE, sep = "\t", comment.char = "#", stringsAsFactors = FALSE)
  
  # Filter for gene entries only
  gene_order <- gtf_data %>%
    filter(V3 == "gene") %>%
    select(V1, V4, V5, V9) %>% # Chromosome, Start, End, Attributes
    filter(str_detect(V9, "gene_id")) %>% # Ensure we only process lines that contain "gene_id"
    mutate(Ensembl_Gene_ID = str_extract(V9, "(ENSG[0-9]+)")) %>% # Use str_extract to capture ENSG number
    select(Ensembl_Gene_ID, V1, V4, V5) %>% # Select relevant columns
    rename(Chromosome = V1, Start = V4, End = V5)
  
  
  # Download chromosome arm information (cytoBand data)
  if (!file.exists(arm_order_file)) {
    download.file(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
                  destfile = arm_order_file)
  }
  
  # Load cytoBand file
  cytoBand <- read.table(gzfile(arm_order_file), header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  colnames(cytoBand) <- c("chrom", "start", "end", "band", "stain")
  
  # Extract chromosome arm information
  cytoBand <- cytoBand %>%
    mutate(arm = substr(band, 1, 1))  # Extract 'p' or 'q' from the band column
  
  # Define arm boundaries
  chromosome_arms <- cytoBand %>%
    group_by(chrom, arm) %>%
    summarize(
      Start = min(start),
      End = max(end),
      .groups = "drop"
    )
  
  # Add chromosome arm information
  gene_order <- gene_order %>%
    # Join chromosome arm information with gene_order based on Chromosome
    left_join(chromosome_arms, by = c("Chromosome" = "chrom")) %>%
    # Determine arm (p or q)
    mutate(Arm = case_when(
      (Start.x > End.y & arm == "p") ~ "q",
      arm == "q" ~ NA,
      .default = arm)) %>%
    na.omit() %>%      
    # Create Chromosome_arm by pasting Chromosome with Arm information
    mutate(Chromosome = paste0(Chromosome, Arm)) %>%
    # Define chromosome arm order
    mutate(Chromosome = factor(Chromosome, levels = c(paste0("chr", rep(1:22, each = 2), c("p", "q")),
                                                      "chrXp", "chrXq", "chrYp", "chrYq"))) %>%
    # Sort genes by Chromosome arm and Start position
    arrange(Chromosome, Start.x)  %>%
    # Select only relevant column for infercnv
    select(Ensembl_Gene_ID, Chromosome, Start.x, End.x) %>%
    # Remove ENSG duplicated (genes that are both on X and Y chromosome need to be remove before infercnv)
    distinct(Ensembl_Gene_ID, .keep_all = TRUE)
  
  # Save the final output
  write_tsv(gene_order, gene_arm_order_file, col_names = FALSE, append = FALSE)

  }
