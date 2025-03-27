This directory contains scripts used in the analysis.

* `00-make-gene-order-file.R` prepares the input gene order file needed for `inferCNV`
* `build-normal-reference/` contains scripts that prepare normal references for use for `inferCNV` for individual projects:
  * `build-reference-SCPCP000015.R` builds references for `SCPCP000015`, including one reference with all immune cells and one with a subset of immune cells (macrophage and T cell types)
