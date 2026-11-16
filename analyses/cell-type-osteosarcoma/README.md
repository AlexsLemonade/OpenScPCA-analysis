# Osteosarcoma cell type annotation

## Description

The goal of this module is to annotate osteosarcoma samples across ScPCA projects, including `SCPCP000017`, `SCPCP000018`, and `SCPCP000023`. 

As described in [this GitHub Discussion](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/1410), this module uses the [`OsteoCAR`](https://figshare.com/articles/dataset/OsteoCAR_A_multi-species_single-cell_atlas_of_primary_and_metastatic_osteosarcoma/31029559) reference to annotate samples by:

* Annotating samples with `SingleR`
* Annotating samples with `scANVI`
* Reconciling these annotations with consensus cell type annotations to derive a final set of annotations across osteosarcoma samples

## Usage

The module can be run with the `run-analysis.sh` script:

```sh
# run full analysis
./run-analysis.sh

# run with test data
testing=1 ./run-analysis.sh
```

## Input files

This module requires the processed `SCE` objects for projects `SCPCP0000{17,18,23}`.

These files can be obtained with the following code run from the root of the `OpenScPCA-analysis` repository.
You must be logged into your AWS account to download these files.

```sh
# Download processed SCE objects
./download-data.py --projects SCPCP000017,SCPCP000018,SCPCP000023 --format SCE
```

## Output files

`run-analysis.sh` creates the following files in the `references` directory:

* `patient_harmonized_aggr_ref.rds`: an aggregated reference for annotating samples with `SingleR`, created from the patient primary and metastasis `OsteoCAR` references.
  Labels are harmonized across both references, using `Ann_Level1` labels for Tumor and Mesenchymal cells and `Ann_Level2` labels for all other cell types.
  When using this reference with `SingleR`, use the `label` column as the reference labels.
**Note:** in the xenograft references, all cells are labeled `Tumor` in `Ann_Level1`, so this reference has a single `Tumor` label and should not be used on its own for cell typing with `SingleR`.

## Software requirements

This module will use both `renv` and `conda` to manage R and Python software environments, respectively.

The `conda` environment is forthcoming. 


To set up the `renv` environment, launch `R` from the module directory and run:

```sh
renv::restore()
```

## Computational resources

TBD
