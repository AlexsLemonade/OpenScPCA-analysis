# Wilms Tumor Dataset Annotation (SCPCP000006) 

Wilms tumor (WT) is the most common pediatric kidney cancer characterized by an exacerbated intra- and inter- tumor heterogeneity. The genetic landscape of WT is very diverse in each of the histological contingents. The COG classifies WT patients into two groups: the favorable histology and diffuse anaplasia. Each of these groups is composed of the blastemal, epithelial, and stromal populations of cancer cells in different proportions, as well as cells from the normal kidney, mostly kidney epithelial cells, endothelial cells, immune cells and normal stromal cells (fibroblast).

## Description

Here, we first aim to annotate the Wilms Tumor snRNA-seq samples in the SCPCP000006 (n=40) dataset. To do so we will:

• Provide annotations of normal cells composing the kidney, including normal kidney epithelium, endothelium, stroma and immune cells
• Provide annotations of tumor cell populations that may be present in the WT samples, including blastemal, epithelial, and stromal populations of cancer cells
Based on the provided annotation, we would like to additionally provide a reference of marker genes for the three cancer cell populations, which is so far lacking for the WT community.

The analysis is/will be divided as the following:

[x] Metadata file: compilation of a metadata file of marker genes for expected cell types that will be used for validation at a later step
[x] Script: clustering of cells across a set of parameters for few samples
[ ] Script: label transfer from the fetal kidney atlas reference using runAzimuth
[ ] Script: run InferCNV
[ ] Notebook: explore results from steps 2 to 4 for about 5 to 10 samples
[ ] Script: compile scripts 2 to 4 in a RMardown file with required adjustements and render it across all samples
[ ] Notebook: explore results from step 6, integrate all samples together and annotate the dataset using (i) metadatafile, (ii) CNV information, (iii) label transfer information

## Usage
From Rstudio, run the Rmd reports or render the R scripts (see below R studio session set up). Please before running the script, make sure that the paths are correct. 
You can also simply have a look at the html reports in the notebook folder. Here, no need to run anything, we try to guide you through the analysis. Have a look at the code using the unhide code button on the top right of each chunk!

## Input files

### Single cell data

In this module, we start with the processed `SingleCellExperiment` objects from the ScPCA Portal.
Data have been downloaded locally and are found in mnt_data. the mnt_data folder has to be define in the config.yaml file or changed in the notebook accordingly. 

```{r paths}
path_to_data <- "~/mnt_data/Wilms ALSF/SCPCP000006_2024-06-25"
```

We choosed to re-build a Seurat object from the counts data and to follow the Seurat workflow [normalization –> reduction –> clustering] 

We transferred meta.data from the _processed.rds object to keep:

- QC data computed by the DataLab

- annotation data computed by the DataLab

- raw annotation and gene_symbol conversion


### metadata 

The SCPCP000006_metadata.tsv file in cell-type-wilms-tumor-06 contains clinical information related to the samples in the dataset. Some information can be helpful for annotation and validation:

- treatment: Some of the samples have been pre-treated with chemotherapy and some are upfront resection. We expect few changes between the 2 conditions, including a higher immune infiltration and more DNA damages pathways in treated samples.

- histology: the COG classifies Wilms tumor as either (i) Favorable or (ii) Anaplastic. Some differenices are expected, some marker genes or pathways are associated with anaplasia (see sets of marker gene). 



## Output files

### Marker sets 

This folder is a resource for later validation of the annotated cell types.

#### The table CellType_metadata.csv contains the following column and information:
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


#### The table GeneticAlterations_metadata.csv contains the following column and information:
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


### 01-clusering

We uploaded a notebook 01-clustering_SCPCS000169.Rmd and html report in the notebook folder. 
This a template for an analysis notebook using R Markdown.
In this notebook, we set up parameters in the Seurat workflow [normalization –> reduction –> clustering] for one Wilms tumor sample (SCPCS000169) of the Wilms Tumor dataset (SCPCP000006) and try to get a first feeling a cells composing the sample. 

After discussing with the DataLab, this template will be adapted and rendered to the 40 Wilms tumor samples. 

## Software requirements

To perform the analysis, run the RMarkdown script in R (version 4.4.1).
The main packages used are:
- Seurat version 5
- Azimuth version 5
- inferCNV
- SCpubr for visualization
- DT for table visualization
- DElegate for differential expression analysis

For complete reproducibility of the results, you can build and run the docker image using the Dockerfile. This will allow you to work on RStudio (R version 4.4.1) from the based image bioconductor/tidyverse:3.19.

In the config.yaml file, define your system specific parameter and paths (e.g. to the data).
Execute the run.sh file and open RStudio in your browser (http://localhost:8080/). 
By default, username = rstudio, password = wordpass.




## Computational resources

No need to run any analysis, just open the metadata table!
