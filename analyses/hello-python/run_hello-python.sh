#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Execute and render the python notebook to html
jupyter nbconvert --execute --to html hello.ipynb
