# Wilms Tumor Dataset Annotation (SCPCP000006)

Wilms tumor (WT) is the most common pediatric kidney cancer characterized by an exacerbated intra- and inter- tumor heterogeneity.
The genetic landscape of WT is very diverse in each of the histological contingents.
The COG classifies WT patients into two groups: the favorable histology and diffuse anaplasia.
Each of these groups is composed of the blastemal, epithelial, and stromal populations of cancer cells in different proportions, as well as cells from the normal kidney, mostly kidney epithelial cells, endothelial cells, immune cells and normal stromal cells (fibroblast).

## Description

Here, we first aim to annotate the Wilms Tumor snRNA-seq samples in the SCPCP000006 (n=40) dataset. To do so we will:

• Provide annotations of normal cells composing the kidney, including normal kidney epithelium, endothelium, stroma and immune cells
• Provide annotations of tumor cell populations that may be present in the WT samples, including blastemal, epithelial, and stromal populations of cancer cells
Based on the provided annotation, we would like to additionally provide a reference of marker genes for the three cancer cell populations, which is so far lacking for the WT community.

The analysis is/will be divided as the following:

- [x] Metadata file: compilation of a metadata file of marker genes for expected cell types that will be used for validation at a later step
- [x] Script: clustering of cells across a set of parameters for few samples
- [x] Script: label transfer from the fetal kidney atlas reference using runAzimuth
- [ ] Script: run InferCNV
- [ ] Notebook: explore results from steps 2 to 4 for about 5 to 10 samples
- [ ] Script: compile scripts 2 to 4 in a RMardown file with required adjustements and render it across all samples
- [ ] Notebook: explore results from step 6, integrate all samples together and annotate the dataset using (i) metadatafile, (ii) CNV information, (iii) label transfer information

## Usage
From Rstudio, run the Rmd reports or render the R scripts (see below R studio session set up).
You can also simply have a look at the html reports in the notebook folder.
Here, no need to run anything, we try to guide you through the analysis. Have a look at the code using the unhide code button on the top right of each chunk!

## Input files

### single nuclei data

We work with the _processed.rds SingleCellExperiment objects.
From the module directory, make sure that the conda environment is set-up:

```shell
conda activate openscpca
```

log into AWS CLI:
```shell
# replace `openscpca` with your AWS CLI profile name if it differs
export AWS_PROFILE=openscpca
aws sso login
```

use download-data.py to download the data as the following:
```shell
../../download-data.py --projects SCPCP000006
```
This is saving the data in OpenScPCA-analysis/data/current/SCPCP000006

Of note, this requires AWS CLI setup to run as intended: https://openscpca.readthedocs.io/en/latest/technical-setup/environment-setup/configure-aws-cli/

### sample metadata

The OpenScPCA-analysis/data/current/SCPCP000006/single_cell_metadata.tsv file contains clinical information related to the samples in the dataset.
Some information can be helpful for annotation and validation:

- treatment: Some of the samples have been pre-treated with chemotherapy and some are upfront resection.
We expect few changes between the 2 conditions, including a higher immune infiltration and more DNA damages pathways in treated samples.

- histology: the COG classifies Wilms tumor as either (i) Favorable or (ii) Anaplastic.
Some differenices are expected, some marker genes or pathways are associated with anaplasia (see sets of marker gene).

## Output files

for each of the steps, we have two types of `output`:

- the `notebook` saved in the `notebook` directory, with a subfolder for each sample. 

- the created objects saved in `results` directory, with a subfolder for each sample. 


# Analysis

## Marker sets

We first build a resource for later validation of the annotated cell types. 
We gather from the litterature marker genes and specific genomic alterations that could help us characterizing the Wilms tumor ecosystem, including cancer and non-cancer cells. 

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
|11p15|loss|malignant|NA|10.1128/mcb.9.4.1799-1803.1989|NA|NA|
|16q|loss|malignant|NA|NA|1317258|Associated_with_relapse|
|1p|loss|malignant|NA|NA|8162576|Associated_with_relapse|
|1q|gain|malignant|NA|10.1016/S0002-9440(10)63982-X|NA|Associated_with_relapse|


## Clustering and label transfer from fetal references

R Script to be rendered : `00_run_workflow.R`

### Introduction

The `00_run_workflow.R` contains the following steps:

- define paths

- download and create the fetal kidney reference: `download-and-create-fetal-kidney-ref.R` in `scripts`

- characterize the fetal kidney reference: `00b_characterize_fetal_kidney_reference_Stewart.Rmd` in `notebook_template`

- loop for each samples:

-- `Seurat workflow`, nornalization and clustering: `01_seurat-processing.Rmd` in `notebook_template`
-- `Azimuth` label transfer from the fetal full reference (Cao et al.) in `notebook_template`
-- `Azimuth` label transfer from the fetal kidney reference (Stewart et al.) in `notebook_template`

### Justification 

The use of the right reference is crucial. 
It is recommended that the cell types in the reference is representative to the cell types to be annotated in the query.

Wilms tumors can contain up to three histologies that resemble fetal kidney: blastema, stroma, and epithelia [1-2].
Because of their histological similarity to fetal kidneys, Wilms tumors are thought to arise from developmental derangements in embryonic renal progenitors.

We thus decided to test and compare two fetal (kidney) references that could be use in the analysis module.

##### Human fetal kidney atlas Stewart et al.

We first wanted to try the human fetal kidney atlas to transfer label into the Wilms tumor samples using azimuth. 
You can find more about the human kidney atlas here: https://www.kidneycellatlas.org/ [3]

##### Human Azimuth fetal reference from Cao et al.

Azimuth also provide a human fetal atlas as a reference [4]. 

The data can be found on Zenodo: 
https://zenodo.org/records/4738021#.YJIW4C2ZNQI

The reference contain cells from 15 organs including kidney from fetal samples. 
Here we will use `Azimuth` to transfer labels from the reference.

### Input and outputs

We start with the `_process.Rds` data to run `01_seurat-processing.Rmd`. 
The output of `01_seurat-processing.Rmd` is saved in `results` in a subfolder for each sample and is the input of the second step `02a_label-transfer_fetal_full_reference_Cao.Rmd`.
The output of `02a_label-transfer_fetal_full_reference_Cao.Rmd` is then the input of `02b_label-transfer_fetal_kidney_reference_Stewart.Rmd`.

At the end of the workflow, we have a `Seurat`object that contains:
- normalization and clustering, dimensional reductions
- label transfer from the fetal full reference
- label transfer from the fetal kidney reference

## Software requirements

To perform the analysis, run the RMarkdown script in R (version 4.4.1).
The main packages used are:
- Seurat version 5
- Azimuth version 5
- inferCNV
- SCpubr for visualization
- DT for table visualization
- DElegate for differential expression analysis

### Docker

To build the Docker image, run the following from this directory:

```shell
docker buildx build . -t openscpca/cell-type-wilms-tumor-06
```

The image will also be available from ECR: <https://gallery.ecr.aws/openscpca/cell-type-wilms-tumor-06>

To run the container and develop in RStudio Server, run the following **from the root of the repository**, Replacing `{PASSWORD}`, including the curly braces, with a password of your choosing:

```shell
docker run \
  --mount type=bind,target=/home/rstudio/OpenScPCA-analysis,source=$PWD \
  -e PASSWORD={PASSWORD} \
  -p 8787:8787 \
  public.ecr.aws/openscpca/cell-type-wilms-tumor-06:latest
```

This will pull the latest version of the image from ECR if you do not yet have a copy locally.

Navigate to <http://localhost:8787/> and log in with the username `rstudio` and the password you set.

Within RStudio Server, `OpenScPCA-analysis` will point to your local copy of the repository.

#### A note on Apple Silicon

If you are on a Mac with an M series chip, you will not be able to use RStudio Server if you are using a `linux/amd64` or `linux/x86_84` (like the ones available from ECR).
You must build an ARM image locally to be able to use RStudio Server within the container.

#### A note for Halbritter lab internal development
This work has been developed on a system that uses podman instead of docker. The steps to run the docker/podman images are slightly different and we saved in run-podman-internal.sh our internal approach to run the container. Please, refer to the Docker section to build and run the container instead. 

### renv

This module uses `renv`.
If you are using RStudio Server within the container, the `renv` project will not be activated by default.
You can install packages within the container and use `renv::snapshot()` to update the lockfile without activating the project without a problem in our testing.
The `renv` lockfile is used to install R packages in the Docker image.

## Computational resources


## References 

- [1] https://www.ncbi.nlm.nih.gov/books/NBK373356/ 

- [2] https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9915828/ 

- [3] https://www.science.org/doi/10.1126/science.aat5031 

- [4] https://www.science.org/doi/10.1126/science.aba7721

