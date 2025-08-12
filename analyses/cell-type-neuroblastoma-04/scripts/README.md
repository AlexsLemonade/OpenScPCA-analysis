This directory contains scripts used in the `cell-type-neuroblastoma-04` module.

* `00_convert-nbatlas.R` converts an [Seurat `NBAtlas` object](https://data.mendeley.com/datasets/yhcf6787yp/3) into `SingleCellExperiment` (`SCE`) and `AnnData` objects for eventual use with `SingleR` and `scANVI`, respectively
* `01_train-singler-model.R` trains a `SingleR` model from an `NBAtlas` reference object, which may or may not be aggregated first
* `02_classify-singler.R` performs cell type annotation with `SingleR` using a given trained `SingleR` model on a given `SCE` object
* `03a_train-scanvi-model.py` trains a `scANVI` model on the `NBAtlas` reference
* `03b_prepare-scanvi-query.R` subsets and prepares a given SCPCP000004 `SCE` object for input to `03c-run-scanvi-label-transfer.py`
* `03c-run-scanvi-label-transfer.py` performs cell type annotation using label transfer using the `NBAtlas` `scANVI` model on a query AnnData object prepared by `03b_prepare-scanvi-query.R`

## utils directory

This directory contains files with utility functions used in code throughout the module.

* `jaccard-utils.R` contains utility functions to build heatmaps colored by Jaccard similarity
* `celltype-utils.R` contains additional utility functions to plot and wrangle cell type results


## setup directory

The `setup` directory contains additional scripts used to prepare for analysis.
These scripts are not included in `../run-analysis.sh`; usage information is provided in scripts directly.

* `prepare-nbatlas-gene-lists.R` is manually run to create `../references/nbatlas-marker-genes.tsv`
