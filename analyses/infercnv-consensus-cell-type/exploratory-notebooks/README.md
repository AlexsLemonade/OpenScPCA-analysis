This directory contains exploratory notebooks that are not specifically run in the analysis pipeline.
Exploratory notebooks for pooled normal reference analyses are organized in directories named for the ScPCA project they analyze.


* `reference-query-counts.Rmd`: This sample looks at the distribution of reference-intended cells to determine how many samples could be run through `inferCNV`
TODO: This notebook will need to be updated once `SCimilarity` annotations are incorporated into consensus cell types.

## `SCPCP000004`

* `01_nb-consensus-cell-types.Rmd` explores the distribution of normal cells across samples in `SCPCP000004` to determine approaches for creating a pooled normal reference


## `SCPCP000015`

* `01_ewings-consensus-cell-types.Rmd` explores the distribution of immune cells across samples in `SCPCP000015` to determine approaches for creating a pooled normal reference
* `02_ewings-reference-cnv.Rmd` explores the distribution of inferred total CNV in pooled normal reference cells across samples in `SCPCP000015`
* `03_ewings-pooled-internal.Rmd` compares results from `SCPCP000015` samples that were able to be processed with both internal and pooled normal references
