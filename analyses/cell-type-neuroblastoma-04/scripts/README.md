This directory contains scripts used in the `cell-type-neuroblastoma-04` module.

* `00a_extract-nbatlas-ids.R` extracts cell ids into a text file from the `seuratObj_NBAtlas_share_v20241203.rds` NBAtlas file for use when converting the reference
* `00b_convert-nbatlas.R` converts an [Seurat `NBAtlas` object](https://data.mendeley.com/datasets/yhcf6787yp/3) into `SingleCellExperiment` and `AnnData` objects for eventual use with `SingleR` and `scANVI`, respectively
