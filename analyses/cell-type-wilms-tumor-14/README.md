# Wilms tumor annotation (SCPCP000014)

We plan to annotate cell types for the Wilms tumor samples (n=10) in [SCPCP000014](https://scpca.alexslemonade.org/projects/SCPCP000014). Our analysis involves data clean-up for low quality nuclei and doublets, cell type annotation, and tumor cell identification.

## Description

The goal of this analysis is to curate a validated cell type annotation for Wilms tumor samples in the portal (SCPCP000014). Specifically, we aim to generate following outcomes: (i) Lists of marker genes to identify cell types in Wilms tumor; (ii) Identification of tumor cells from normal cells; (iii) Refined annotation of cell types among normal cells; (iv) Annotation of sub-groups among tumor cells, if applicable.

#### 00. Pre-processing the provided SCE objects
This would include:
* formatting SCE objects to Seurat objects for further analysis
* count normalization, feature selection, transformation, PCA, UMAP, batch effect correction (if merged object)

#### 01. Anchor transfer using Seurat
  * Fetal reference from [Human Kidney atlas](https://www.kidneycellatlas.org/) works best in my preliminary analysis.

#### 02. Curating marker gene lists
- Tumor cell markers for Wilms tumor
- Kidney cell types

#### 03. Cell type annotation with marker gene lists
* Cellassign
* scType

#### 04. Tumor cell identification
- inferCNV?
- CopyKat?
- Based on Tumor marker genes

#### 05. Sample merging and validation

## Usage

* Run main pipeline with command line
```bash
cd /path/to/OpenScPCA-analysis
cd analyses/cell-type-wilms-tumor-14
bash ./run_cell-type-wilms-14.sh
```

## Input files

* Download processed SCE objects
```bash
# download wilms tumor dataset (n=10)
cd /path/to/OpenScPCA-analysis
./download-data.py --projects SCPCP000014
./download-results.py --modules doublet-detection --projects SCPCP000014
```

## Output files

All results are sync under S3 bucket `researcher-009160072044-us-east-2`.

#### 00. Pre-processing the provided SCE objects 
- Results for this section contains 10 `.rdsSeurat` objects for further analysis, stored in `scratch/` folder.

#### 01. Anchor transfer using Seurat
- Path on S3 `s3://researcher-009160072044-us-east-2/cell-type-wilms-tumor-14/results/01_anchor_transfer_seurat`. 
- For more description, see `./results/README.md`.

## Software requirements

- Install missing dependencies on AWS virtual computer:
```bash
sudo apt install -y libglpk40 \
  libcurl4-openssl-dev \
  jags \
  libmagick++-dev \
  libhdf5-dev \
  libharfbuzz-dev libfribidi-dev libtiff5-dev
```
- Versions for required R packages listed in `./renv.lock`

## Computational resources

Analysis could be executed on a virtual computer ([Standard-4XL](https://openscpca.readthedocs.io/en/latest/aws/lsfr/creating-vcs/)) via AWS Lightsail for Research.

## Exploratory analysis
In addition to the main pipeline, some exploratory analysis in R notebooks are added into the `./exploratory_analysis` folder, including CNV analysis. Check `./exploratory_analysis/README.md` for more details.