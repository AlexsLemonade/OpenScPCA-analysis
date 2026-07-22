# Osteosarcoma cell type annotation

## Description

The goal of this module is to annotate osteosarcoma samples across ScPCA projects, including `SCPCP000017`, `SCPCP000018`, and `SCPCP000023`. 

As described in [this GitHub Discussion](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/1410), this module uses the [`OsteoCAR`](https://figshare.com/articles/dataset/OsteoCAR_A_multi-species_single-cell_atlas_of_primary_and_metastatic_osteosarcoma/31029559) reference to annotate samples by:

* Annotating samples with `SingleR`
* Annotating samples with `scANVI`
* Reconciling these annotations with consensus cell type annotations to derive a final set of annotations across osteosarcoma samples

## Usage

## Input files

This module requires the processed `SCE` objects for projects `SCPCP0000{17,18,23}`.

These files can be obtained with the following code run from the root of the `OpenScPCA-analysis` repository.
You must be logged into your AWS account to download these files.

```sh
# Download processed SCE objects
./download-data.py --projects SCPCP000017,SCPCP000018,SCPCP000023 --format SCE
```

## Output files

TBD

## Software requirements

This module will both `renv` and `conda` to manage R and Python software environments, respectively.

Actual environments are forthcoming. 

## Computational resources

TBD
