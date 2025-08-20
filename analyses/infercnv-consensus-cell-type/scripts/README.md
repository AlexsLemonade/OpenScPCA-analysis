This directory contains scripts used in the analysis.

* `prepare-diagnosis-map.R` is manually run to prepare the TSV file mapping broad diagnosis groups to specific ScPCA diagnoses
* `00-make-gene-order-file.R` prepares the input gene order file needed for `inferCNV`
* `01_run-infercnv.R` runs `inferCNV` on a given ScPCA library using a specified normal reference
  * The normal reference can either be "pooled" (created from combining cells across project libraries) or "internal" (created from cells in the given library only)
* `utils.R` contains helper functions used by other scripts in this directory

## Normal references

The script `build-normal-references.R` can be used to build a normal reference for a given project.
Note that this was developed after processing `SCPCP000015`, whose references were built with `build-reference-SCPCP000015.R`.

Normal reference files should be stored in `../references/normal-references/{project id}` and named as `{reference cell group}.rds`.
These SCE objects should contain at least:

* A raw counts matrix
* Column names (cell ids) formatted as `{library_id}-{barcode}`
* A `colData` column `consensus_annotation` recording each cell's consensus annotation
