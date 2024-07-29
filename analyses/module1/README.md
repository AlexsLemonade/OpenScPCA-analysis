# Module 1: Compilation of a metadata file of marker genes for expected cell types

Wilms tumor (WT) is the most common pediatric kidney cancer characterized by an exacerbated intra- and inter- tumor heterogeneity. The genetic landscape of WT is very diverse in each of the histological contingents. The COG classifies WT patients into two groups: the favorable histology and diffuse anaplasia. Each of these groups is composed of the blastemal, epithelial, and stromal populations of cancer cells in different proportions, as well as cells from the normal kidney, mostly kidney epithelial cells, endothelial cells, immune cells and normal stromal cells (fibroblast).

## Description

In this module, we reviewed the literature to compile a table of marker genes for each of the expected cell types in the dataset. 

## Usage

This module is a resource for later validation of the annotated cell types. The table of marker genes can be found in the folder results under the name CellType_metadata.csv

## Output files

The table contains the following column and information:
- "gene_symbol" contains the symbol of the described gene, using the HUGO Gene Nomenclature
- ENSEMBL_ID contains the stable identifier from the ENSEMBL database
- cell_class is either "malignant" for marker genes specific to malignant population, or "non-malignant" for markers genes specific to non-malignant tissue or "both" for marker genes that can be found in malignant as well as non-malignant tissue but are still informative in respect to the cell type.
- cell_type contains the list of the cell types that are attributed to the marker gene
- DOI contains the list of main publication identifiers supporting the choice of the marker gene
- comment can be empty or contains any additional information


## Software requirements

No software required.

## Computational resources

No need to run any analysis, just open the metadata table!
