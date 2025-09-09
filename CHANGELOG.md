# Release notes for OpenScPCA-analysis

This document contains release notes for versioned releases of the OpenScPCA-analysis repository.
We will also periodically release and tag dated snapshots of the repository, but these will not be reflected in this document.

<!--
Add new release notes in reverse numerical order (newest first) below this comment

You may want to add temporary notes here for tracking as features are added, before a new release is ready.
-->

## v0.2.3

The following includes the following new modules developed by the Data Lab:  

- `infercnv-consensus-cell-type`: a module for exploring different strategies to use consensus cell types as a normal reference when running [`InferCNV`](https://github.com/broadinstitute/inferCNV/wiki) 
- `cell-type-neuroblastoma-04`: a module for annotating the Neuroblastoma samples in `SCPCP000004`
- `cell-type-scimilarity`: a module for annotating cells using [`SCimilarity`](https://genentech.github.io/scimilarity/index.html)

The `cell-type-consensus` module was updated to now incorporate annotations from the `cell-type-scimilarity` module such that consensus cell types are now assigned if 2 out of 3 methods agree (`SingleR`, `CellAssign`, and `SCimilarity`). 

This release includes minor updates to the following modules: 

- `cell-type-wilms-06` 
- `cell-type-ETP-ALL-03`
- `cell-type-nonETP-ALL-03`
- `cell-type-ewings`

The following community-contributed analysis modules have been initiated and are in active development: 

- `cell-type-osteosarcoma-23`

## v0.2.2

This release includes minor updates to the following modules: 

- `cell-type-consensus`
- `cell-type-ewings`
- `cell-type-wilms-06`

Within the `cell-type-consensus` module, the [reference used to assign consensus cell types](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/analyses/cell-type-consensus/references/consensus-cell-type-reference.tsv) has been updated to now include `neuron` as a possible cell type. 

## v0.2.1

This release includes minor updates to the following modules: 

- `cell-type-ewings`
- `cell-type-consensus`

In particular, the workflow for assigning consensus cell types has been created in `cell-type-consensus` module. 
This workflow is now present in `OpenScPCA-nf`, and results will be made available to all contributors. 

New FAQs covering how to work with ScPCA data as `Seurat` objects and using gene symbols instead of Ensembl IDs have been added to the [OpenScPCA documentation](https://openscpca.readthedocs.io/en/latest/troubleshooting-faq/faq/). 

`Rcpp` package versions within modules were updated for compatibility. 

## v0.2.0

This release adds the first set of community-contributed analyses to the repository.
These modules are focused on cell type identification and annotation for specific ScPCA datasets.
Note that many of these modules are still in development at this stage and may not be fully functional.

- `cell-type-dsrct`
- `cell-type-ETP-ALL-03`
- `cell-type-glioblastoma`
- `cell-type-nonETP-ALL-03`
- `cell-type-wilms-tumor-06`
- `cell-type-wilms-tumor-14`


This release also adds the following new modules developed by the Data Lab:

- `cell-type-consensus`: a module for exploring consensus cell types across multiple annotation methods
- `hello-clusters`: a demonstration module for clustering analysis using the [`rOpenScPCA` package](https://github.com/AlexsLemonade/rOpenScPCA)
- `seurat-conversion`: a module for converting `SingleCellExperiment` objects to Seurat objects, also using the `rOpenScPCA` package

Other updates in this release include:

- a new `sync-results.py` script to simplify uploading (and downloading) analysis results from an analysis module to a user's S3 bucket
- changes from `miniconda` to `miniforge` for conda usage throughout the project

While not part of this repository, we do want to also note that we have created the [`rOpenScPCA` package](https://github.com/AlexsLemonade/rOpenScPCA), which will house utility functions commonly used by analysis modules here.
The goal is to centralize common functions used across analysis modules to make it easier to share code and maintain consistency across modules.


## v0.1.0

This is the initial versioned release of the OpenScPCA-analysis repository.
The repository at this stage should be generally complete with respect to infrastructure, but the analysis is still only in the early stages.

With respect to infrastructure, the repository contains the following components:

- detailed documentation in the `docs` directory explaining how to interact with the OpenScPCA project and how to set up and run analyses (published at https://openscpca.readthedocs.io)
- a `create-analysis-module.py` script for setting up new analysis modules
- `download-data.py` and `download-results.py` scripts to download data and results from the OpenScPCA project
- template notebooks, scripts, environment files, and Docker images for analysis modules

The repository contains the following analysis modules that are now in a relatively mature state:

- `hello-R`
- `hello-python`
- `simulate-sce`
- `doublet-detection`

The following additional modules are in active development:

- `cell-type-ewings`
- `metacells`
