# Marker genes

This folder contains a table of marker genes for the cell types in the aggregated OsteoCAR reference, along with the script used to add Ensembl gene IDs to the table.

## Files

* `osteosarcoma-markers.tsv`: marker genes for each cell type, with one row per gene and the following columns:
  * `cell_type`: the cell type label, matching the labels in the aggregated reference
  * `gene_symbol`: the gene symbol
  * `ensembl_gene_id`: the Ensembl gene ID
  * `source`: the DOI of the reference where the gene is reported for that cell type
* `convert-Ensembl-ID.R`: adds Ensembl gene IDs to `osteosarcoma-markers.tsv` based on the gene symbols, using `org.Hs.eg.db`

## Usage

To add or update the Ensembl gene IDs, run the script from the module directory:

```sh
Rscript marker-genes/convert-Ensembl-ID.R
```

The script prints any gene symbols that match no Ensembl ID or more than one, so those can be checked manually.
