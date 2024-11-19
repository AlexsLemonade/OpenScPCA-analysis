This directory contains supplementary notebooks that are not used as part of this module's workflow.

- `compare-label-transfer-approaches.Rmd` compares label transfer results for 5 samples from two approaches: i) using Azimuth directly, and ii) using code adapted from Azimuth functions.
  - This notebook was written because Azimuth was not able to be used with OpenScPCA test data (as described in <https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/810>).
  We therefore adapted code from Azimuth to perform label transfer directly.
  This notebook compares label transfer results from this new approach to those obtained with Azimuth to ensure the labels are reasonably similar.