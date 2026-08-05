#!/bin/bash

# This script runs the osteosarcoma cell type annotation module.
#
# Usage:
#
# ./run-analysis.sh
#
# There are several variables that can be defined when calling this script:
# - testing (Default value: 0)
#   - Use `testing=1` to run with test data
#   - Example usage: testing=1 ./run-analysis.sh
# - `threads` (Default value: 4)
#   - Use `threads=X` to specify a different number of threads
#   - Example usage: threads=2 ./run-analysis.sh # to request 2 threads

set -euo pipefail

# Ensure script is being run from its directory
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}

# Define and create directories
script_dir="scripts"
ref_dir="references"
mkdir -p $ref_dir

# Define argument defaults
testing=${testing:-0} # default is not testing
threads=${threads:-4} # default 4 threads


# Download OsteoCar reference ----------------------------------------------
figshare_mets_url="https://api.figshare.com/v2/file/download/65376111"
figshare_prim_url="https://api.figshare.com/v2/file/download/65375988"
figshare_xeno_mets_url="https://api.figshare.com/v2/file/download/65375973"
figshare_xeno_prim_url="https://api.figshare.com/v2/file/download/65375976"
mets_qs2="${ref_dir}/patient_mets.qs2"
prim_qs2="${ref_dir}/patient_prim.qs2"
xeno_mets_qs2="${ref_dir}/xeno_mets.qs2"
xeno_prim_qs2="${ref_dir}/xeno_prim.qs2"

# First, download the object with a helper function
# This function takes two arguments in order, the URL and the filename to save to
download_file() {
  local url="$1"
  local file_name="$2"
  if [[ ! -f $file_name ]]; then
    curl -Lo $file_name $url
  fi
}

# Download OsteoCar reference files
download_file $figshare_mets_url $mets_qs2
download_file $figshare_prim_url $prim_qs2
download_file $figshare_xeno_mets_url $xeno_mets_qs2
download_file $figshare_xeno_prim_url $xeno_prim_qs2


# Convert OsteoCar to SCE (TODO: and later, also AnnData)
mets_sce="${ref_dir}/patient_mets_sce.rds"
prim_sce="${ref_dir}/patient_prim_sce.rds"
xeno_mets_sce="${ref_dir}/xeno_mets_sce.rds"
xeno_prim_sce="${ref_dir}/xeno_prim_sce.rds"


Rscript ${script_dir}/convert-osteocar.R --input_ref_file $mets_qs2 --output_sce_file $mets_sce
Rscript ${script_dir}/convert-osteocar.R --input_ref_file $prim_qs2 --output_sce_file $prim_sce
Rscript ${script_dir}/convert-osteocar.R --input_ref_file $xeno_mets_qs2 --output_sce_file $xeno_mets_sce
Rscript ${script_dir}/convert-osteocar.R --input_ref_file $xeno_prim_qs2 --output_sce_file $xeno_prim_sce
