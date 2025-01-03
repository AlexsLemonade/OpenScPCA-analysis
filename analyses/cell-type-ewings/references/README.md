# Reference files

This directory contains reference files used for cell typing Ewing sarcoma samples from SCPCP000015.

## Tumor marker genes

`tumor-marker-genes.tsv` contains a list of marker genes used to identify tumor cells in previously published studies.
Currently, this list only contains genes that should mark all tumor cells and we are not yet considering subsets of tumor cells.

The marker gene list was created using gene lists from two publications:

1. Visser, Bleijs, _et al._ (2023) https://doi.org/10.1158/2767-9764.CRC-23-0027 - Marker genes were obtained from the list in [Supplemental Methods under `Cell annotation using external data`](https://aacr.silverchair-cdn.com/aacr/content_public/journal/cancerrescommun/3/10/10.1158_2767-9764.crc-23-0027/1/crc-23-0027-s01.pdf).

2. Goodspeed _et al._ (2024) https://doi.org/10.1101/2024.01.18.576251 - Marker genes were determined to be those that were unique to Ewing sarcoma cells in Figure 1C and 1F.
Only genes that were consistently expressed across Ewing sarcoma tumor cells with little expression in non-malignant cell types were chosen as marker genes.

## All marker genes

`visser-all-marker-genes.tsv` contains a list of marker genes used to identify tumor cells and normal cells.
This currently contains all genes in `tumor-marker-genes.tsv` and genes used to identify normal cells in [Visser, Beligs, _et al._ (2023)](https://doi.org/10.1158/2767-9764.CRC-23-0027) as indicated in the supplemental methods.

## Panglao marker genes

`PanglaoDB_markers_2020-03-27.tsv` is the full set of marker genes exported from [PanglaoDB](https://panglaodb.se), versioned at the given date

## CellAssign References

The `cellassign_refs` folder contains any binary matrices that are used as references when running `CellAssign`.

These files were created in `01-marker-gene-classification-cellassign.Rmd`:

1. `tumor-marker_cellassign.tsv`: This file contains a binary matrix with all genes in `tumor-marker-genes.tsv`.
2. `filtered-tumor-marker_cellassign.tsv`: This file contains a binary matrix with only the genes that have mean gene expression > 1 in `SCPCS000490`.
3. `visser-all-marker_cellassign.tsv`: This file contains a binary matrix with all marker genes in `tumor-marker-genes.tsv` and all markers for normal cells identified in the Supplemental methods of [Visser et al.,](https://doi.org/10.1158/2767-9764.CRC-23-0027).
4. `panglao-endo-fibro_cellassign.tsv`: This file contains a binary matrix with all marker genes for endothelial cells and fibroblasts from `PanglaoDB_markers_2020-03-27.tsv` and all tumor markers from `tumor-marker-genes.tsv`.

## InferCNV references

The `infercnv_refs` folder contains any references needed to run `InferCNV`.

`InferCNV` requires a gene order file containing all genes and the start and stop positions for those genes.
This file is a tab delimited `.txt` file with no column headers.
The columns correspond to Ensembl gene id, chromosome, start, and stop.

The gene order file, `Homo_sapiens.GRCh38.104.gene_order.txt`, is too large to store on the remote repository.
If you need to run any analysis for `InferCNV`, this file can be created by running `scripts/make-gene-order-file.R` with default parameters.

The other required input for `InferCNV` is an annotation file with two columns and no column headers.
The first column contains the cell barcode and the second contains the annotation for that cell (either `reference` or `unknown`).
These files are specific for each library and depend on which cells are denoted as the reference.
Each library contains a folder with any annotations file used to run `InferCNV` for that library.

## Marker gene sets for identifying tumor cell states 

The `tumor-cell-state-markers.tsv` file contains a list of marker genes that can be used to classify tumor cell states in Ewing samples. 
The marker genes included here are specific to EWS-FLI1 high, EWS-FLI1 low, and proliferative tumor cells. 
This list was obtained based on key genes mentioned in [Goodspeed _et al._](https://doi.org/10.1101/2024.01.18.576251), [Aynaud _et al._](https://doi.org/10.1016/j.celrep.2020.01.049), and [Wrenn _et al._](https://doi.org/10.1158/1078-0432.CCR-23-1111). 

### Gene signatures 

The `gene_signatures` folder contains any custom gene lists obtained from publications that can be used to identify tumor cell states: 

1. `anyaud-ews-targets.tsv`: A list of the 78 marker genes defined by [Aynaud _et al._](https://doi.org/10.1016/j.celrep.2020.01.049) to be EWS-FLI1 targets. 
Figure 4 shows that expression of these targets is correlated with EWS-FLI1 levels at a single-cell level. 

2. `wrenn-nt5e-genes.tsv`: A list of 28 genes from [Wrenn _et al._](https://doi.org/10.1158/1078-0432.CCR-23-1111) that represent the overlap between the top 217 genes correlated with _NT5E_ expression in patient tumors and the top 200 markers of _NT5E+_ Ewing sarcoma cells _in vitro_. 
These genes are shown in Figure 5D and 5E. 

The following gene sets from MsigDb were also used to define EWS-FLI1 targets and may be helpful in defining cell states:  

- [STAEGE_EWING_FAMILY_TUMOR](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/STAEGE_EWING_FAMILY_TUMOR.html)
- [MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_UP](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_UP.html)
- [MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_DN](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/MIYAGAWA_TARGETS_OF_EWSR1_ETS_FUSIONS_DN.html)
- [ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/ZHANG_TARGETS_OF_EWSR1_FLI1_FUSION.html)
- [RIGGI_EWING_SARCOMA_PROGENITOR_UP](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/RIGGI_EWING_SARCOMA_PROGENITOR_UP.html) 
- [RIGGI_EWING_SARCOMA_PROGENITOR_DN](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/RIGGI_EWING_SARCOMA_PROGENITOR_DN.html) 
- [KINSEY_TARGETS_OF_EWSR1_FLI1_FUSION_UP](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_UP.html)
- [KINSEY_TARGETS_OF_EWSR1_FLI1_FUSION_DN](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/KINSEY_TARGETS_OF_EWSR1_FLII_FUSION_DN.html)

Wrenn _et al._ also used found that the following additional gene sets were highly expressed in CD73 high, EWS-FLI1 low tumor cells: 

- [HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.html)
GO:BP ECM Organization 
- [GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION.html)
