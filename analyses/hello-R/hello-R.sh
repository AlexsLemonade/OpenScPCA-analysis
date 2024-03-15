#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Render the R notebook

notebook="notebooks/hello.Rmd"

Rscript -e "rmarkdown::render('${notebook}')"
