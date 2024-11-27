This directory contains supplementary notebooks that are not as part of this module's workflow.

`characterize_fetal_kidney_reference_Stewart.Rmd` explores the Stewart et al. fetal kidney reference.
  - This notebook exports results from this exploration to `../results/reference/`, but these are not directly used in the workflow.
- `compare-label-transfer-approaches.Rmd` compares label transfer results for 5 samples from two approaches: i) using Azimuth directly, and ii) using code adapted from Azimuth functions.
  - This notebook was written because Azimuth was not able to be used with OpenScPCA test data (as described in <https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/810>).
  We therefore adapted code from Azimuth to perform label transfer directly.
  This notebook compares label transfer results from this new approach to those obtained with Azimuth to ensure the labels are reasonably similar.
- `cnv-exploration`
  - This directory contains exploratory notebooks that are only run on a subset of samples to explore CNV methods.
  - These notebooks run generated with the script `../scripts/explore-cnv-methods.R`.
  - Rendered HTML notebooks are stored in sample-specific subfolders, and they are not fully up-to-date with their R Markdown counterparts and/or the OpenScPCA data release.
