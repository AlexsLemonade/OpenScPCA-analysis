#!/bin/bash
set -euo pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Render the python notebook

jupyter execute hello.ipynb
