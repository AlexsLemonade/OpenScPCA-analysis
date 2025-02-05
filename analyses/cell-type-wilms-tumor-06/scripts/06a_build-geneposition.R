# This script downloads and adapt the genome position file for use with infercnv

# load library
library(tidyverse)


# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")
# the paths to the results folder in which output references will be saved
reference_dir <- file.path(module_base, "results", "references")

# Define the file paths
gene_order_file <- file.path(reference_dir, "gencode_v38_gene_pos.txt")
arm_order_file <- file.path(reference_dir, "hg38_cytoBand.txt")
gene_arm_order_file <- file.path(reference_dir, "gencode_v38_gene_pos_arm.txt")

if (!file.exists(gene_arm_order_file)) {
  
  # Download the updated gene order file for GENCODE v38 (Ensembl 104)
  # Read the GTF file and extract gene positions
  gtf_data <- read_tsv(gene_order_file, col_names = FALSE, comment = "#")
  colnames(gtf_data) <- c("Chromosome", "HAVANA", "Type", "Start", "End", "X6", "Direction", "X8", "Description")
  
  # Filter for gene entries only
  gene_order <- gtf_data %>%
    filter(Type == "gene") %>%
    filter(str_detect(Description, "gene_id")) %>% # Ensure we only process lines that contain "gene_id"
    mutate(Ensembl_Gene_ID = str_extract(Description, "(ENSG[0-9]+)")) %>% # Use str_extract to capture ENSG number
    select(Ensembl_Gene_ID, Chromosome, Start, End)  # Select relevant columns

  
  # Download chromosome arm information (cytoBand data)
  if (!file.exists(arm_order_file)) {
    download.file(url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz",
                  destfile = arm_order_file)
  }
  
  # Load cytoBand file
  cytoBand <- read_tsv(arm_order_file, col_names = FALSE)
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
    ) %>%
    na.omit() %>%  # remove NA values
    filter(arm == "p")  # we only keep one `arm` here, as the information of interest is the position of the centromere, which is the `End` of the p arm and the `Start` of the q arm
 
   # Add chromosome arm information
  gene_order <- gene_order %>%
    # Join chromosome arm information with gene_order based on Chromosome
    left_join(chromosome_arms, by = c("Chromosome" = "chrom")) %>%
    # Determine arm (p or q)
    # by default and construction of the `chromosome_arms`, arm = `p`
    # if the gene is after the centromere, we change arm to `q`
    mutate(arm = case_when(
      (Start.x > End.y ) ~ "q",
      .default = arm)) %>%
    # Create Chromosome_arm by pasting Chromosome with Arm information
    mutate(Chromosome = paste0(Chromosome, arm))  %>%
    # Define chromosome arm order
    mutate(Chromosome = factor(Chromosome, levels = c(paste0("chr", rep(1:22, each = 2), c("p", "q")),
                                                      "chrXp", "chrXq", "chrYp", "chrYq", "chrM"))) %>%
    # Sort genes by Chromosome arm and Start position
    arrange(Chromosome, Start.x)  %>%
    # Select only relevant column for infercnv
    select(Ensembl_Gene_ID, Chromosome, Start.x, End.x) %>%
    # Remove ENSG duplicated (genes that are both on X and Y chromosome need to be remove before infercnv)
    distinct(Ensembl_Gene_ID, .keep_all = TRUE)
  
  # Save the final output
  write_tsv(gene_order, gene_arm_order_file, col_names = FALSE, append = FALSE)

  }
