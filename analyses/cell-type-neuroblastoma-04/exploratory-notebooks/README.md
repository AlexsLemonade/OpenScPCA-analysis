This directory contains exploratory notebooks that are not specifically run in the analysis pipeline.

* `reference-aggregation.Rmd` is a template notebook that performs a cursory comparison of SingleR inferences made with an aggregated vs non-aggregated reference for a given library, as well as a brief comparison to consensus cell types
  * Rendered notebooks are stored in `reference-aggregation_htmls` and were generated with the script `generate_reference-aggregation.sh`
* `separate-tumor.Rmd` is a template notebook that performs a cursory comparison of SingleR inferences made with an (aggregated) reference that considers either a single `Neuroendocrine` category or two categories `Neuroendocrine` and `Neuroendocrine-tumor` where the latter contains cells in the NBAtlas "tumor zoom"
  * Rendered notebooks are stored in `separate-tumor_htmls` and were generated with the script `generate_separate-tumor.sh`
* `jaccard-utils.R` contains utility functions to build heatmaps colored by Jaccard similarity
