This directory contains exploratory notebooks which are not run by `../run-analysis.sh`.

* `01-check-SCPCP000023.Rmd` compares samples in `SCPCP000023` to the `OsteoCar` reference to determine which samples are represented in both
* `02-compare-references.Rmd` compares profiles from `OsteoCar` references to one another to begin understanding their contents and relationship, for the goal of identifying which to use
* `03a-compare-patient-references-ann-level2.Rmd` compares some cursory results across several types of reference constructions by running `SingleR` on two patient samples which also appear in `OsteoCar`, focusing on the `Ann_Level2` labels
  * This notebook uses functions provided in `exploratory-utils.Rmd`
* `03b-compare-patient-references-combined-ann-levels.Rmd` performs a similar analysis as `03a-compare-patient-references-ann-level2.Rmd`, except using the `Ann_Level1` donations for Tumor and Mesenchymal cell types but the `Ann_Level2` labels for the other categories
  * This notebook uses functions provided in `exploratory-utils.Rmd`
 
