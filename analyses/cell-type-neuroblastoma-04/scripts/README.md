This directory contains scripts used in the `cell-type-neuroblastoma-04` module.

* `00a_extract-nbatlas-ids.R` extracts cell ids into a text file from the `seuratObj_NBAtlas_share_v20241203.rds` NBAtlas file for use when converting the reference
* `00b_convert-nbatlas.R` converts an [Seurat `NBAtlas` object](https://data.mendeley.com/datasets/yhcf6787yp/3) into `SingleCellExperiment` and `AnnData` objects with Ensembl gene ids for eventual use with `SingleR` and `scANVI`, respectively
  * Note that some genes are necessarily lost from the converted objects where there is no equivalent ScPCA identifier, or where duplicate gene symbols cannot be resolved
