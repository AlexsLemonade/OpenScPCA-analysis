This directory contains scripts used in the `cell-type-neuroblastoma-04` module.

* `00_convert-nbatlas.R` converts an [Seurat `NBAtlas` object](https://data.mendeley.com/datasets/yhcf6787yp/3) into `SingleCellExperiment` and `AnnData` objects for eventual use with `SingleR` and `scANVI`, respectively
* `01_train-singler-model.R` trains a SingleR model from an `NBAtlas` reference object, which may or may not be aggregated first
* `02_classify-singler.R` performs cell type annotation with SingleR using a given trained SingleR model on a given SCE object

## Setup directory

The `setup` directory contains additional scripts used to prepare for analysis.
These scripts are not included in `../run-analysis.sh`; usage information is provided in scripts directly.

* `prepare-gene-lists.R` is manually run to create `../references/nbatlas-marker-genes.tsv` and `../references/consensus-marker-genes.tsv` validation files
