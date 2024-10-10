# we download and adapt the genome position file

# The base path for the OpenScPCA repository, found by its (hidden) .git directory
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
# The path to this module
module_base <- file.path(repository_base, "analyses", "cell-type-wilms-tumor-06")

# load library
library(tidyverse)
# Define the gene order file
  gene_order_file <- file.path(module_base, "results","references", 'gencode_v19_gene_pos.txt')

# Download the file if ot existing
if (!file.exists(gene_order_file)) {
  download.file(url = 'https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gen_pos.complete.txt',
                destfile = gene_order_file)
  gene_order <- read_tsv(gene_order_file, col_names =  FALSE)
  tmp <- gsub("\\..*","",gene_order$X1) 
  tmp <- gsub(".*\\|","",tmp) 
  gene_order$X1 <- tmp
  gene_order <- gene_order[grepl("ENSG", x = gene_order$X1),]
  
  write_tsv(col_names = FALSE, gene_order, gene_order_file, append = FALSE)
}
  