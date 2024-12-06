#!/bin/bash

# This script renders this module's example notebooks.
# This script is primarily used for testing.
# Please refer to the individual notebooks in this module for examples of performing and evaluating clustering with rOpenScPCA.

set -euo pipefail

# Ensure script is being run from its directory
MODULE_DIR=$(dirname "${BASH_SOURCE[0]}")
cd ${MODULE_DIR}

# Render first notebook
Rscript -e "rmarkdown::render('01_perform-evaluate-clustering.Rmd', clean = TRUE)"
