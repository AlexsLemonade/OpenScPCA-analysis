#!/bin/bash

# This script is used to download the cell marker reference file 
# The downloaded file will be saved to `references/Cell_marker_Human.xlsx`

cell_marker_file="${scripts_dir}/../references/Cell_marker_Human.xlsx"
cell_marker_url="http://117.50.127.228/CellMarker/CellMarker_download_files/file/Cell_marker_Human.xlsx"

if [[ ! -f $cell_marker_file ]]; then 
    curl -o $cell_marker_file $cell_marker_url
fi
