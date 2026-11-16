# Read file containing marker genes for osteosarcoma
marker_file <- "marker-genes/osteosarcoma-markers.tsv"

# Read the existing marker table, preserving its source citations.
markers <- readr::read_tsv(marker_file, show_col_types = FALSE)

# map each gene symbol to its Ensembl gene ID(s)
ensembl_ids <- AnnotationDbi::mapIds(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = markers$gene_symbol,
  keytype = "SYMBOL",
  column = "ENSEMBL",
  multiVals = "list"
)

# check for symbols with no ID or with more than one ID
ensembl_ids[lengths(ensembl_ids) != 1]

# use the first ID for each gene (after checking any with more than one)
markers$ensembl_gene_id <- purrr::map_chr(ensembl_ids, \(ids) ids[1])

# save the final file
readr::write_tsv(markers, "marker-genes/osteosarcoma-markers.tsv")
