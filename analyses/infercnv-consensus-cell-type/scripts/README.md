This directory contains scripts used in the analysis.

* `00-make-gene-order-file.R` prepares the input gene order file needed for `inferCNV`
* `01_run-infer-cnv.R` runs `inferCNV` on a given ScPCA library


## Normal references

The directory `build-normal-reference/` contains scripts that prepare normal SCE references for use for `inferCNV` for individual projects and are named `build-reference-{project id}.R`.
Normal reference files should be stored in `../references/normal-references/` and named as `ref-{reference name}.rds`.
These SCE objects should contain at least:

* A raw counts matrix
* Column names (cell ids) formatted as `{library_id}-{barcode}`
* A `colData` column `consensus_annotation` recording each cell's consensus annotation
