# Wilms tumor annotation (SCPCP000014)

We plan to annotate cell types for the Wilms tumor samples (n=10) in [SCPCP000014](https://scpca.alexslemonade.org/projects/SCPCP000014). Our analysis involves data clean-up for low quality nuclei and doublets, cell type annotation, and tumor cell identification.

## Description

The goal of this analysis is to curate a validated cell type annotation for Wilms tumor samples in the portal (SCPCP000014). Specifically, we aim to generate following outcomes: (i) Lists of marker genes to identify cell types in Wilms tumor; (ii) Identification of tumor cells from normal cells; (iii) Refined annotation of cell types among normal cells; (iv) Annotation of sub-groups among tumor cells, if applicable.

#### 00. Pre-processing the provided SCE objects
This would include:
* formatting SCE objects to Seurat objects for further analysis
* count normalization, feature selection, transformation, PCA, UMAP, batch effect correction (if merged object)

#### 01. Anchor transfer using Seurat

#### 02. Curating marker gene lists
- Tumor cell markers for Wilms tumor
- Kidney cell types

#### 03. Cell type annotation with marker gene lists
* Cellassign
* scType

#### 04. Tumor cell identification
- inferCNV (no reference, confused by results)
- CopyCat?
- Based on Tumor marker genes

#### 05. Sample merging and validation

## Usage

* Run scripts interactively on Rstudio.
```R
Rscript --vanilla xx.R
```
Please provide instructions on how to run the analysis module.
What commands are needed to execute all steps in the analysis?

## Input files

* Download processed SCE objects
```bash
# download wilms tumor dataset (n=10)
cd /path/to/OpenScPCA-analysis
./download-data.py --projects SCPCP000014
./download-results.py --modules doublet-detection --projects SCPCP000014
```
* Create this module structure
```bash
./create-analysis-module.py cell-type-wilms-SCPCP000014 --use-r  --use-renv --use-conda --conda-file-only
```

## Output files

Please include a description of the output from your analysis, including:

- What type of files are created?
- What are the contents of the output files?
- Where are the files stored?
- Are any intermediate files generated?
If so, where are they stored?

## Software requirements

- Setup conda channel priority
```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

- Create conda env and other packages not available on conda
```bash
cd /home/lightsail-user/git/OpenScPCA-analysis/analyses/cell-type-wilms-tumor-14
conda env create -f ./conda_envs/main.yml -y -n wilms-tumor-14-main
conda activate wilms-tumor-14-main
Rscript --vanilla ./conda_envs/main.R
```

## Computational resources

Analysis could be executed on a virtual computer ([Standard-2XL](https://openscpca.readthedocs.io/en/latest/aws/lsfr/creating-vcs/)) via AWS Lightsail for Research.