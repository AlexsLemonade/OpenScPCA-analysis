# Module 1: Compilation of a metadata file of marker genes for expected cell types

Wilms tumor (WT) is the most common pediatric kidney cancer characterized by an exacerbated intra- and inter- tumor heterogeneity. The genetic landscape of WT is very diverse in each of the histological contingents. The COG classifies WT patients into two groups: the favorable histology and diffuse anaplasia. Each of these groups is composed of the blastemal, epithelial, and stromal populations of cancer cells in different proportions, as well as cells from the normal kidney, mostly kidney epithelial cells, endothelial cells, immune cells and normal stromal cells (fibroblast).

## Description

In this module, we reviewed the literature to compile a table of marker genes for each of the expected cell types in the dataset. Additionally, we provide a table of know genetic alterations in Wilms tumor that can be useful to validate CNV profiles obtained after running inferCNV function. 

## Usage

This module is a resource for later validation of the annotated cell types. The table of marker genes can be found in the folder results under the name CellType_metadata.csv. The table of known genetic alterations can also be found in the results folder under the name GeneticAlterations_metadata.csv.

## Output files

### The table CellType_metadata.csv contains the following column and information:
- "gene_symbol" contains the symbol of the described gene, using the HUGO Gene Nomenclature
- ENSEMBL_ID contains the stable identifier from the ENSEMBL database
- cell_class is either "malignant" for marker genes specific to malignant population, or "non-malignant" for markers genes specific to non-malignant tissue or "both" for marker genes that can be found in malignant as well as non-malignant tissue but are still informative in respect to the cell type.
- cell_type contains the list of the cell types that are attributed to the marker gene
- DOI contains the list of main publication identifiers supporting the choice of the marker gene
- comment can be empty or contains any additional information

  |gene_symbol|ENSEMBL_ID|cell_class|cell_type|DOI|comment|
  |---|---|---|---|---|---|
  |WT1|ENSG00000184937|malignant|cancer_cell|10.1242/dev.153163|Tumor_suppressor_WT1_is_lost_in_some_WT_cells|
  |IGF2|ENSG00000167244|malignant|cancer_cell|10.1038/ng1293-408|NA|
  |TP53|ENSG00000141510|malignant|anaplastic|10.1158/1078-0432.CCR-16-0985|Might_also_be_in_small_non_anaplastic_subset|
  |MYCN|ENSG00000134323|malignant|anaplastic|10.18632/oncotarget.3377|Also_in_non_anaplastic_poor_outcome|
  |MAX|ENSG00000125952|malignant|anaplastic|10.1016/j.ccell.2015.01.002|Also_in_non_anaplastic_poor_outcome|
  |SIX1|ENSG00000126778|malignant|blastema|10.1016/j.ccell.2015.01.002|NA|
  |SIX2|ENSG00000170577|malignant|blastema|10.1016/j.ccell.2015.01.002|NA|
  |CITED1|ENSG00000125931|malignant|blastema|10.1593/neo.07358|Also_in_embryonic_kidney|
  |PTPRC|ENSG00000081237|immune|NA|10.1101/gr.273300.120|NA|
  |CD68|ENSG00000129226|immune|myeloid|10.1186/1746-1596-7-12|NA|
  |CD163|ENSG00000177575|immune|macrophage|10.1186/1746-1596-7-12|NA|
  |VWF|ENSG00000110799|endothelium|endothelium|10.1134/S1990747819030140|NA|
  |CD3E|ENSG00000198851|immune|T_cell|10.1101/gr.273300.120|NA|
  |MS4A1|ENSG00000156738|immune|B_cell|10.1101/gr.273300.120|NA|
  |FOXP3|ENSG00000049768|immune|T_cell|10.1101/gr.273300.120|Treg|
  |CD4|ENSG00000010610|immune|T_cell|10.1101/gr.273300.120|NA|
  |CD8A|ENSG00000153563|immune|T_cell|10.1101/gr.273300.120|NA|
  |EPCAM|ENSG00000119888|NA|epithelial|10.1016/j.stemcr.2014.05.013|epithelial_malignant_and_non_malignant|
  |NCAM1|ENSG00000149294|malignant|blastema|10.1016/j.stemcr.2014.05.013|might_also_be_expressed_in_non_malignant|
  |PODXL|ENSG00000128567|non-malignant|podocyte|10.1016/j.stem.2019.06.009|NA|
  |COL6A3|ENSG00000163359|malignant|mesenchymal|10.2147/OTT.S256654|might_also_be_expressed_in_non_malignant_stroma|
  |THY1|ENSG00000154096|malignant|mesenchymal|10.1093/hmg/ddq042|might_also_be_expressed_in_non_malignant_stroma|





### The table GeneticAlterations_metadata.csv contains the following column and information:
- alteration contains the number and portion of the affected chromosome
- gain_loss contains the information regarding the gain or loss of the corresponding genetic alteration
- cell_class is "malignant" 
- cell_type contains the list of the malignant cell types that are attributed to the marker gene, either blastemal, stromal, epithelial or NA if none of the three histology is more prone to the described genetic alteration
- DOI contains the list of main publication identifiers supporting the choice of the genetic alteration
- comment can be empty or contains any additional information

|alteration|gain_loss|cell_class|cell_type|DOI|PMID|comment
|---|---|---|---|---|---|---|
|11p13|loss|malignant|NA|10.1242/dev.153163|NA|NA|
|11p15|loss|malignant|NA|10.1128/mcb.9.4.1799|NA|NA|
|16q|loss|malignant|NA|NA|1317258|Associated_with_relapse|
|1p|loss|malignant|NA|NA|8162576|Associated_with_relapse|
|1q|gain|malignant|NA|10.1016/S0002-9440(10)63982-X|NA|Associated_with_relapse|



## Software requirements

No software required.

## Computational resources

No need to run any analysis, just open the metadata table!
