#!/bin/bash

# This script runs the analysis module

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

echo "Running 00-01_processing_rds.R"
Rscript scripts/00-01_processing_rds.R

echo "Running 02-03_annotation.R"
Rscript scripts/02-03_annotation.R

echo "Running 04_multipanel_plot.R"
Rscript scripts/04_multipanel_plot.R

echo "Running 05_cluster_evaluation.R"
Rscript scripts/05_cluster_evaluation.R

echo "Running 06_sctype_exploration.R"
Rscript scripts/06_sctype_exploration.R

echo "Running 07_run_copykat.R"
Rscript scripts/07_run_copykat.R

echo "Running markerGenes_submission.R"
Rscript scripts/markerGenes_submission.R

echo "Running writeout_submission.R"
Rscript scripts/writeout_submission.R
