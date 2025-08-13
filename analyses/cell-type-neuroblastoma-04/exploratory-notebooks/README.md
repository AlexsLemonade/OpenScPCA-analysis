This directory contains exploratory notebooks that are not specifically run in the analysis pipeline.


* `singler-results.Rmd` explores `SingleR` annotation results and compares them to consensus cell types
* `scanvi-singler-comparison.Rmd` compares `scANVI` and `SingleR` annotation results to one another, including a brief comparison to consensus cell types


## `singler-testing`

This directory contains notebooks and scripts from exploratory analyses planning how to run `SingleR`.

* `reference-aggregation.Rmd` is a template notebook that performs a cursory comparison of SingleR inferences made with an aggregated vs non-aggregated reference for a given library, as well as a brief comparison to consensus cell types
  * Rendered notebooks are stored in `reference-aggregation_htmls` and were generated with the script `generate_reference-aggregation.sh`
* `separate-tumor.Rmd` is a template notebook that performs a cursory comparison of SingleR inferences made with an (aggregated) reference that considers either a single `Neuroendocrine` category or two categories `Neuroendocrine` and `Neuroendocrine-tumor` where the latter contains cells in the NBAtlas "tumor zoom"
  * Rendered notebooks are stored in `separate-tumor_htmls` and were generated with the script `generate_separate-tumor.sh`
* `filter-genes.Rmd` is a template notebook that performs a cursory comparison of SingleR inferences made with a reference where mitochondrial and ribosomal genes were filtered out of the reference vs. not filtered
  * Rendered notebooks are stored in `filter-genes_htmls` and were generated with the script `generate_filter-genes.sh`
