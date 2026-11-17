# Scripts

This folder contains scripts used in this analysis module. 

* `00-download-metafactory-results.sh`: This script is used to download and save results from running `metafactory` on ScPCA samples. 
Running this script will save results to the `scratch` directory in this module to be used as input to other scripts and notebooks. 
Currently, only Data Lab members with appropriate AWS access can run this. 
Be sure to set your `AWS_PROFILE` environment variable before running. 

* `neuroblastoma-exploratory-plots.R`: This script was used to create two figures used in the poster and presentation on `metafactory` presented at the 2026 Pediatric AACR meeting.
