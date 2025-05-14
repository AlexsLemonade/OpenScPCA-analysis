This directory contains scripts used in the analysis.

* `00-make-gene-order-file.R` prepares the input gene order file needed for `inferCNV`
* `01_run-infercnv.R` runs `inferCNV` on a given ScPCA library using a specified normal reference
  * The normal reference can either be "pooled" (created from combining cells across project libraries) or "internal" (created from cells in the given library only)
* `utils.R` contains helper functions used by other scripts in this directory

## Normal references

The directory `build-normal-references/` contains scripts that prepare pooled normal SCE references for use with `inferCNV` for individual projects and are named `build-reference-{project id}.R`.
Normal reference files should be stored in `../references/normal-references/` and named as `{reference name}.rds`.
These SCE objects should contain at least:

* A raw counts matrix
* Column names (cell ids) formatted as `{library_id}-{barcode}`
* A `colData` column `consensus_annotation` recording each cell's consensus annotation
