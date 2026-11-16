# Read file containing marker genes for osteosarcoma
marker_file <- "marker-genes/osteosarcoma-markers.tsv"

# Make sure the file has the correct format (tabs)
# read the lines, skipping the header, and split each one on any whitespace
lines <- readr::read_lines(marker_file)[-1]
parts <- strsplit(trimws(lines), "\\s+")

markers <- data.frame(
  cell_type = purrr::map_chr(parts, 1),
  gene_symbol = purrr::map_chr(parts, 2),
  ensembl_gene_id = NA,
  source = purrr::map_chr(parts, 3)
)

# check before saving
head(markers)

# save the file with the updated format
readr::write_tsv(markers, marker_file, na = "")

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
