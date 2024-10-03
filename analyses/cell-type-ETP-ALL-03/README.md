# ETP T-ALL Annotation (SCPCP000003)

This analysis module will include codes to annotate cell types and tumor/normal status in ETP T-ALL from SCPCP000003 (n=30) present on the ScPCA portal.

## Description

We first aim to annotate the cell types in ETP T-ALL, and use the annotated B cells in the sample as the "normal" cells to identify tumor cells, since T-ALL is caused by the clonal proliferation of immature T-cell [<https://www.nature.com/articles/s41375-018-0127-8>].

-   We use the cell type marker (`Azimuth_BM_level1.csv`) from [Azimuth Human Bone Marrow reference](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Bone%20Marrow). In total, there are 14 cell types: B, CD4T, CD8T, Other T, DC, Monocytes, Macrophages, NK, Early Erythrocytes, Late Erythrocytes, Plasma, Platelet, Stromal, and Hematopoietic Stem and Progenitor Cells (HSPC). Based on the exploratory analysis, we believe that most of the cells in these samples do not express adequate markers to be distinguished at finer cell type level (eg. naive vs memory, CD14 vs CD16 etc.), and majority of the cells should belong to T-cells. In addition, we include the marker genes for blast cell [[Bhasin et al. (2023)](https://www.nature.com/articles/s41598-023-39152-z)] as well as erythroid precursor and cancer cell in immune system [[ScType](https://sctype.app/database.php) database].

-   Since ScType annotates cell types at cluster level using marker genes provided by user or from the built-in database, we employ [self-assembling manifold](https://github.com/atarashansky/self-assembling-manifold/tree/master) (SAM) algorithm, a soft feature selection strategy for better separation of homogeneous cell types.

-   After cell type annotation, we provide B cells as the normal cells in the sample, if there is any, to [CopyKat](https://github.com/navinlabcode/copykat), for identification of tumor cells.

Here are the steps in the module:

1.  Generating a processed rds file for each sample using SAM (`scripts/00-01_processing_rds.R`)

2.  Annotating cell type using ScType and identifying tumor cells using CopyKat (`scripts/02-03_annotation.R`)

## Usage

Before running Rscripts in R or Rstudio, we first need to prepare the input files as shown in the next section, and run the following codes in the terminal for installing required libraries:

```         
#system packages installation
sudo apt install libglpk40
sudo apt install libcurl4-openssl-dev                 #for Seurat
sudo apt-get install libxml2-dev libfontconfig1-dev libharfbuzz-dev  libfribidi-dev libtiff5-dev  #for devtools

conda-lock install --name openscpca-cell-type-ETP-ALL-03 conda-lock.yml
Rscript -e "renv::restore()"
```

## Input files

The `scripts/00-01_processing_rds.R` requires the processed SingleCellExperiment objects (`_processed.rds`) and doublet-detection results (`_processed_scdblfinder.tsv`) from SCPCP000003. These files could be obtained from running the following codes:

```         
#run in terminal
../../download-data.py --projects SCPCP000003
../../download-results.py --projects SCPCP000003 --modules doublet-detection
```

As for the annotation, `scripts/02-03_annotation.R` requires cell type marker gene file, `Azimuth_BM_level1.csv`, as an input for ScType. This excel file contains a list of positive marker genes in Ensembl ID under `ensembl_id_positive_marker` for each cell type; *TMEM56* and *CD235a* are not detected in our dataset, thus they are being removed as part of the markers for Late Eryth and Pre Eryth respectively. As of now, there is no negative marker genes provided under `ensembl_id_negative_marker`.

## Output files

Running `scripts/00-01_processing_rds.R` will generate two types of output:

-   `rds` objects in `scratch/`

-   umap plots showing leiden clustering in `plots/`

Running `scripts/02-03_annotation.R` will generate several outputs:

-   updated `rds` objects in `scratch/`

-   umap plots showing cell type and CopyKat prediction (if there is any) and dotplots showing the features added with `AddModuleScore()` in `plots/`

-   ScType results of top 10 possible cell types in a cluster (`_sctype_top10_celltypes_perCluster.txt`) and metadata file tabulating leiden cluster, cell type, low confidence cell type, and CopyKat prediction for each cell (`_metadata.txt`) in `results/`

## Software requirements

To run the analysis, execute the Rscript in R or Rstudio (version 4.4.0). The main libraries used are:

-   Seurat (version 5.1.0)

-   reticulate (version 1.39.0)

-   sam-algorithm (in python)

-   ScType

-   CopyKat

The renv.lock file contains all packages and version information. All python libraries are installed in the conda environment `openscpca-cell-type-ETP-ALL-03`, and the python codes are executed in the same environment by running them in R via `reticulate`. To create and activate this environment from `.yml` file use:

```         
conda-lock install --name openscpca-cell-type-ETP-ALL-03 conda-lock.yml
```

## Computational resources

All the commands above are currently executed in the standard 4XL virtual machine via AWS Lightsail for Research, but it runs pretty slow for CopyKat with one computational core.
