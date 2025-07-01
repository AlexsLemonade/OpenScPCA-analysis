This directory contains scripts used in the `cell-type-neuroblastoma-04` module.

* `00a_extract-nbatlas-ids.R` extracts cell ids into a text file from the `seuratObj_NBAtlas_share_v20241203.rds` NBAtlas file for use when converting the reference
* `00b_convert-nbatlas.R` converts an [Seurat `NBAtlas` object](https://data.mendeley.com/datasets/yhcf6787yp/3) into `SingleCellExperiment` and `AnnData` objects for eventual use with `SingleR` and `scANVI`, respectively
* `01_train-singler-model.R` aggregates given `NBAtlas` reference object and trains an associated SingleR model

## Setup directory

The `setup` directory contains additional scripts used to prepare for analysis.
These scripts are not included in `../run-analysis.sh`; usage information is provided in scripts directly.

* `prepare-gene-lists.R` is manually run to create `../references/nbatlas-marker-genes.tsv` and `../references/consensus-marker-genes.tsv` validation files
