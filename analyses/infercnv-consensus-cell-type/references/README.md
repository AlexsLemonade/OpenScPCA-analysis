This directory contains reference information needed for analysis.

* `cnv-validation.tsv`: Lists specific CNVs which can be used for validation for diagnoses explored in the module
* `Homo_sapiens.GRCh38.104.gene_order.txt` gene order file for input to `inferCNV`
  * This file is created with `../scripts/00-make-gene-order-file.R` when running the analysis pipeline, but it is ignored from the repository due to its large size
* `normal-references/`: This directory contains, organized by project, RDS files containing SCE objects for use as normal references with `inferCNV`.
More information about these files is provided below.
This directory is ignored from the repository due to its large size

## Normal references

Normal references for use with `inferCNV` are formatted as SCE files and are expected to have the following components:

* Files should be named as `ref_{reference name}.rds`
* Reference SCEs should contain at least:
  * A raw counts matrix
  * Column names (cell ids) formatted as `{library_id}-{barcode}`
  * A `colData` column `consensus_annotation` recording each cell's consensus annotation


The `normal-references/` directory contains the following references for each given project:

* `SCPCP000015`
  * `ref_immune.rds` contains all immune cells across samples, excluding those which were annotated as `tumor` in the `cell-type-ewings` module
  * `ref_endo.rds` contains all endothelial cells across samples, excluding those which were annotated as `tumor` in the `cell-type-ewings` module
  * `ref_endo-immune.rds` contains the union of `ref_immune` and `ref_endo`