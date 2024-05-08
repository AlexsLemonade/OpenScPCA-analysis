# Reference files

This directory contains reference files used for cell typing Ewing sarcoma samples from SCPCP000015.
All files saved here can also be found in `s3://researcher-211125375652-us-east-2/ewings-cell-type/references`.

Note that some files may be too big to be stored in the remote repository, and are only in the S3 results bucket.
All references will need to be synced to the local repository before using any scripts that require larger references.

## Tumor marker genes

`tumor-marker-genes.tsv` contains a list of marker genes used to identify tumor cells in previously published studies.
Currently, this list only contains genes that should mark all tumor cells and we are not yet considering subsets of tumor cells.

The marker gene list was created using gene lists from two publications:

1. Visser, Bleijs, _et al._ (2023) https://doi.org/10.1158/2767-9764.CRC-23-0027 - Marker genes were obtained from the list in [Supplemental Methods under `Cell annotation using external data`](https://aacr.silverchair-cdn.com/aacr/content_public/journal/cancerrescommun/3/10/10.1158_2767-9764.crc-23-0027/1/crc-23-0027-s01.pdf).

2. Goodspeed _et al._ (2024) https://doi.org/10.1101/2024.01.18.576251 - Marker genes were determined to be those that were unique to Ewing sarcoma cells in Figure 1C and 1F.
Only genes that were consistently expressed across Ewing sarcoma tumor cells with little expression in non-malignant cell types were chosen as marker genes.

## All marker genes

`all-marker-genes.tsv` contains a list of marker genes used to identify tumor cells and normal cells.
This currently contains all genes in `tumor-marker-genes.tsv` and genes used to identify normal cells in [Visser, Beligs, _et al._ (2023)](https://doi.org/10.1158/2767-9764.CRC-23-0027) as indicated in the supplemental methods.

## Panglao marker genes

`PanglaoDB_markers_2020-03-27.tsv` is the full set of marker genes exported from [PanglaoDB](https://panglaodb.se), versioned at the given date

## CellAssign References

The `cellassign_refs` folder contains any binary matrices that are used as references when running `CellAssign`.

Both of these files were created in `01-marker-gene-classification-cellassign.Rmd`:

1. `tumor-marker-cellassign.tsv`: This file contains a binary matrix with all genes in `tumor-marker-genes.tsv`.
2. `filtered-tumor-marker-cellassign.tsv`: This file contains a binary matrix with only the genes that have mean gene expression > 1 in SCPCS000490.
3. `visser-all-marker-cellassign.tsv`: This file contains a binary matrix with all marker genes in `tumor-marker-genes.tsv` and all markers for normal cells identified in the Supplemental methods of [Visser et al.,](https://doi.org/10.1158/2767-9764.CRC-23-0027).
4. `panglao-endo-fibro-cellassign.tsv`: This file contains a binary matrix with all marker genes for endothelial cells and fibroblasts from `PanglaoDB_markers_2020-03-27.tsv` and all tumor markes from `tumor-marker-genes.tsv`.

## InferCNV references

The `infercnv_refs` folder contains any references needed to run `InferCNV`.

`InferCNV` requires a gene order file contianing all genes and the start and stop positions for those genes.
This file is a tab delimited `.txt` file with no column headers.
The columns correspond to Ensembl gene id, chromosome, start, and stop.

The gene order file, `Homo_sapiens.GRCh38.104.gene_order.txt`, was created using the `scripts/make-gene-order-file.R` script with
the reference gtf file, `Homo_sapiens.GRCh38.104.gtf.gz`, as input.
**Note that these files are not present in the remote repository and must be synced using `scripts/sync-references.R` prior to running scripts that use these references.**

The other required input for `InferCNV` is an annotation file with two columns and no column headers.
The first column contains the cell barcode and the second contains the annotation for that cell (either `reference` or `unknown`).
These files are specific for each library and depend on which cells are denoted as the reference.
Each library contains a folder with any annotations file used to run `InferCNV` for that library.
