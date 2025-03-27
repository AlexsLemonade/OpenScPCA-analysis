This directory contains reference information needed for analysis.

* `cnv-validation.tsv`: Lists specific CNVs which can be used for validation for diagnoses explored in the module
* `Homo_sapiens.GRCh38.104.gene_order.txt` gene order file for input to `inferCNV`
  * This file is created with `../scripts/00-make-gene-order-file.R` when running the analysis pipeline, but it is ignored from the repository due to its large size
* `normal-references/`: This directory contains, organized by project, RDS files containing SCE objects for use as normal references with `inferCNV`.
  It is ignored from the repository due to its large size

