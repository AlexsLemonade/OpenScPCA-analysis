# ETP T-ALL Annotation (SCPCP000003)

This analysis module will include codes to annotate cell types and tumor/normal status in ETP T-ALL from SCPCP000003 (n=31) present on the ScPCA portal.

## Description

We first aim to annotate the cell types in ETP T-ALL, and use the annotated B cells in the sample as the "normal" cells to identify tumor cells, since T-ALL is caused by the clonal proliferation of immature T-cell [<https://www.nature.com/articles/s41375-018-0127-8>].

-   We use the cell type marker (`Azimuth_BM_level1.csv`) from [Azimuth Human Bone Marrow reference](https://azimuth.hubmapconsortium.org/references/#Human%20-%20Bone%20Marrow). In total, there are 14 cell types: B, CD4T, CD8T, Other T, DC, Monocytes, Macrophages, NK, Early Erythrocytes, Late Erythrocytes, Plasma, Platelet, Stromal, and Hematopoietic Stem and Progenitor Cells (`HSPC`). Based on the exploratory analysis, we believe that most of the cells in these samples do not express adequate markers to be distinguished at finer cell type level (e.g. naive vs memory, CD14 vs CD16 etc.), and majority of the cells should belong to T-cells. In addition, we include the marker genes for blast cell [[Bhasin et al. (2023)](https://www.nature.com/articles/s41598-023-39152-z)] as well as erythroid precursor and cancer cell in immune system [[`ScType`](https://sctype.app/database.php) database].

    \*\*`Azimuth_BM_level1.csv` is converted to `submission_markerGenes.tsv`, in the final submission format.

-   Since `ScType` annotates cell types at cluster level using marker genes provided by user or from the built-in database, we employ [self-assembling manifold](https://github.com/atarashansky/self-assembling-manifold/tree/master) (SAM) algorithm, a soft feature selection strategy for better separation of homogeneous cell types.

-   After cell type annotation, we fine-tune the annotated B cells by applying 99 percentile cutoff of non-B `ScType` score on the "B cell clusters". We then use the new B cells (i.e those cells which passed the cutoff) as the normal cells in running [`CopyKAT`](https://github.com/navinlabcode/copykat), for the identification of tumor cells.

Here are the steps in the module:

1.  Generating a processed `rds` file for each sample using SAM (`scripts/00-01_processing_rds.R`)

2.  Annotating cell type using `ScType` and identifying tumor cells using `CopyKAT` (`scripts/02-03_annotation.R`)

3.  Fine-tuning the B cells (`scripts/06_sctype_exploration.R`)

4.  Re-running `CopyKAT` (`scripts/07_run_copykat.R`)

## Usage

Before running R scripts in R or RStudio, we first need to prepare the input files as shown in the next section, and run the following codes in the terminal for installing required libraries:

```
#system packages installation
sudo apt install libglpk40
sudo apt install libcurl4-openssl-dev                 #for Seurat
sudo apt-get install libxml2-dev libfontconfig1-dev libharfbuzz-dev  libfribidi-dev libtiff5-dev  #for devtools
sudo apt-get install r-cran-rjags     #for InferCNV, if wish to run

conda-lock install --name openscpca-cell-type-ETP-ALL-03 conda-lock.yml
Rscript -e "renv::restore()"
```

To run the module, use the provided bash script:

```sh
bash run-module.sh
```

## Input files

The `scripts/00-01_processing_rds.R` requires the processed `SingleCellExperiment` objects (`_processed.rds`) and doublet-detection results (`_processed_scdblfinder.tsv`) from SCPCP000003. These files could be obtained from running the following codes:

```
#run in terminal
../../download-data.py --projects SCPCP000003
../../download-results.py --projects SCPCP000003 --modules doublet-detection
```

As for the annotation, `scripts/02-03_annotation.R` requires cell type marker gene file, `Azimuth_BM_level1.csv`, as an input for `ScType`. This excel file contains a list of positive marker genes in Ensembl ID under `ensembl_id_positive_marker` for each cell type; `TMEM56` and `CD235a` are not detected in our dataset, thus they are being removed as part of the markers for `Late Eryth` and `Pre Eryth` respectively. As of now, there is no negative marker genes provided under `ensembl_id_negative_marker`.

## Important output files

-   `rds` objects in `results/rds`

-   `ScType` results of top 10 possible cell types in a cluster (`results/_sctype_top10_celltypes_perCluster.txt`) and `ScType` score (`results/_sctype_scores.txt`)

-   location of fine-tuned B cells in UMAP (`plots/sctype_exploration/_newBcells.png`) and the cell type assignment with added fine-tuned B cells (`results/_newB-normal-annotation.txt`)

-   final submission table (`results/submission_table/_metadata.tsv`) and the UMAP plots showing cell_type_assignment from `ScType` and tumor_cell_classification from `CopyKAT` using fine-tuned B cells (`results/submission_table/multipanels_.png`)

## Software requirements

To run the analysis, execute the R script in R or RStudio (version 4.4.0). The main libraries used are:

-   Seurat (version 5.1.0)

-   `reticulate` (version 1.39.0)

-   `sam-algorithm` (in python)

-   `ScType`

-   `CopyKAT`

The renv.lock file contains all packages and version information. All python libraries are installed in the conda environment `openscpca-cell-type-ETP-ALL-03`, and the python codes are executed in the same environment by running them in R via `reticulate`. To create and activate this environment from `.yml` file use:

```
conda-lock install --name openscpca-cell-type-ETP-ALL-03 conda-lock.yml
```

A Docker image is also available for this module here: <https://gallery.ecr.aws/openscpca/cell-type-etp-all-03>.

## Computational resources

All the commands above are currently executed in the standard 4XL virtual machine via AWS Lightsail for Research, but it runs pretty slow for `CopyKAT` with one computational core.
This module has not been tested on other operating systems.