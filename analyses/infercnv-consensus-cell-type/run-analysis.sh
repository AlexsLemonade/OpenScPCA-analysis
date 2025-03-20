#!/bin/bash

# This script runs the module workflow.
#
# Usage:
#
# ./run-analysis.sh

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

##### Analysis for SCPCP000015 #####
