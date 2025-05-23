# This script builds a TSV of genomic coordinates for input to inferCNV

library(tidyverse)
library(optparse)


# Parse arguments --------------------------------------------------------------
# set up arguments
option_list <- list(
  make_option(
    opt_str = c("--testing"),
    action = "store_true",
    default = FALSE,
    help = "This flag should be specified when test data is being used. The exported file will not include separate rows for chromosome arms."
  )
)
opts <- parse_args(OptionParser(option_list = option_list))

# Setup ------------------

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")
# The path to the reference directory
reference_dir <- file.path(module_base, "results", "references")

# URLs for data to download
gtf_file_url <- "https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz"
arm_order_file_url <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz"


# Define the gene order file
gtf_file <- file.path(reference_dir, "Homo_sapiens.GRCh38.104.gtf.gz")

# Define the arm order file
arm_order_file <- file.path(reference_dir, "hg38_cytoBand.txt")

# Define the gene/arm combined information file for output
gene_arm_order_file <- file.path(reference_dir, "gencode_v38_gene_pos_arm.txt")

# Download files if they don't exist --------------------
if (!file.exists(gtf_file)) {
  download.file(url = gtf_file_url, destfile = gtf_file)
}
if (!file.exists(arm_order_file)) {
  download.file(url = arm_order_file_url, destfile = arm_order_file)
}

# Read and prepare input files -----------------

gtf <- rtracklayer::readGFF(gtf_file)

gene_order_df <- gtf %>%
  # keep only genes on standard chromosomes
  filter(
    grepl("^[0-9XY]+$", seqid),
    type == "gene"
  ) %>%
  select(
    ensembl_id = gene_id,
    chrom = seqid,
    gene_start = start,
    gene_end = end
  ) |>
  mutate(chrom = glue::glue("chr{chrom}"))

# if testing, we export gene_order_df without arm information
if (opts$testing) {
  write_tsv(gene_order_df, gene_arm_order_file, col_names = FALSE, append = FALSE)
} else {
  # Load cytoBand file into R and assign column names
  cytoBand <- read_tsv(arm_order_file, col_names = FALSE)
  colnames(cytoBand) <- c("chrom", "chrom_arm_start", "chrom_arm_end", "band", "stain")

  # Add a column for the chromosome arm (p or q) and group by chromosome and arm
  #  to calculate arm start and end positions
  chromosome_arms_df <- cytoBand %>%
    mutate(arm = substr(band, 1, 1)) %>%
    group_by(chrom, arm) %>%
    summarize(
      chrom_arm_start = min(chrom_arm_start),
      chrom_arm_end = max(chrom_arm_end),
      .groups = "drop"
    ) %>%
    # this will remove non-standard chromosomes which have NA arms
    tidyr::drop_na()

  # Before continuing, we should have 48 rows in chromosome_arms:
  stopifnot("Could not get all chromosome arm bounds" = nrow(chromosome_arms_df) == 48)


  # Combine data frames to get the chromosome and arm for each gene --------------
  combined_df <- gene_order_df %>%
    # combine gene coordinates with chromosome arm coordinates
    left_join(
      chromosome_arms_df,
      by = "chrom",
      relationship = "many-to-many"
    ) %>%
    # keep only rows where gene is on the chromosome arm
    filter(gene_start >= chrom_arm_start & gene_end <= chrom_arm_end) %>%
    # create chrom_arm column
    mutate(chrom_arm = glue::glue("{chrom}{arm}")) %>%
    # Define chromosome arm order
    mutate(chrom_arm = factor(chrom_arm, levels = c(
      paste0("chr", rep(1:22, each = 2), c("p", "q")),
      "chrXp", "chrXq", "chrYp", "chrYq"
    ))) %>%
    # Sort genes by Chromosome arm and Start position
    arrange(chrom_arm, gene_start) %>%
    # Select only relevant column for infercnv
    select(ensembl_id, chrom_arm, gene_start, gene_end)


  # Export --------------
  write_tsv(combined_df, gene_arm_order_file, col_names = FALSE, append = FALSE)
}
