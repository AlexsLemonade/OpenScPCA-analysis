#!/bin/bash

# This script is used to download the reference file for PanglaoDB from scpca-nf
# The downloaded file will be saved to `references/panglao-cell-type-ontologies.tsv`

set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")
scripts_dir=$(pwd)

# define path to ref file and url
ref_file="${scripts_dir}/../references/PanglaoDB_markers_2020-03-27.tsv"
ref_url="https://raw.githubusercontent.com/AlexsLemonade/scpca-nf/refs/heads/main/references/PanglaoDB_markers_2020-03-27.tsv"

# if ref file doesn't exist download from scpca-nf repository 
if [[ ! -f $ref_file ]]; then 
    curl -o $ref_file $ref_url
fi
