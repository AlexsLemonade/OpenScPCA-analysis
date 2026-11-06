#!/bin/bash

# This script is used to download the results from metafactory on ScPCA samples to the scratch directory 
# Currently, only Data Lab staff members have access to these results with the appropriate AWS profile

set -euo pipefail

# navigate to where script lives
cd $(dirname "$0")

# AWS profile 
profile="${AWS_PROFILE:-default}"

# project to use 
project_id="SCPCP000001"

# define inputs and outputs
s3_bucket="s3://ccdl-scpca-nf-results-testing/metafactory/scpca-testing-output/results"

scratch_dir="../scratch"
project_dir="$scratch_dir/$project_id"

# copy results files to local directory 
aws s3 sync "$s3_bucket/$project_id" $project_dir --profile $profile

