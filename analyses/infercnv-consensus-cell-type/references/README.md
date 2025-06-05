This directory contains reference information needed for analysis.

* `cnv-validation.tsv`: TSV file listing specific CNVs which can be used for validation for diagnoses explored in the module
* `Homo_sapiens.GRCh38.104.gene_order.txt`: The gene order file for input to `inferCNV`
  * This file is created with `../scripts/00-make-gene-order-file.R` when running the analysis pipeline, but it is ignored from the repository due to its large size
* `reference-cell-groups.tsv`: TSV file mapping consensus cell types into reference cell group categories to assist in creating normal references across projects
  * The exception to this is `SCPCP000015`, which was analyzed before this TSV file was created
* `normal-references/`: Holds project-specific RDS files containing SCE objects for use as normal references with `inferCNV`.
More information about these files is provided below

## Normal references

Normal references for use with `inferCNV` are formatted as SCE files and are expected to have the following components:

* A raw counts matrix
* Column names (cell ids) formatted as `{library_id}-{barcode}`
* A `colData` column `consensus_annotation` recording each cell's consensus annotation

The `normal-references/` directory contains the following references for each given project:

* `SCPCP000015`
  * `immune.rds` contains all immune cells across non-PDX samples, excluding those which were annotated as `tumor` in the `cell-type-ewings` module
  * `endo.rds` contains all endothelial cells across non-PDX samples, excluding those which were annotated as `tumor` in the `cell-type-ewings` module
  * `endo-immune.rds` contains the union of the `immune` and `endo` references
  * `reference-celltypes.tsv` is a TSV of all consensus cell type labels that contribute to each normal reference
