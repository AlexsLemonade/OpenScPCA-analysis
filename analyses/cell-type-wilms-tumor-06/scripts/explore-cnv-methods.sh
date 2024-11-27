#!/bin/bash

# This script performs exploratory analysis by running copyKAT and inferCNV
# on several samples under different conditions, including with and without a
# specified reference, and exploring their results in notebooks.
#
# Default usage:
# ./explore-cnv-methods.sh
#
# By default, copyKAT will use 32 threads.
# This can be overridden with the THREADS variable:
# THREADS=16 ./explore-cnv-methods.sh
#
# If running with test data, set the variable
# TESTING=1 to ensure that inferCNV is run
# without a reference, since test data will not
# reliably have normal cells to use for the reference.
# TESTING=1 ./explore-cnv-methods.sh
#

set -euo pipefail

THREADS=${THREADS:-32}
TESTING=${TESTING:-0}

# Ensure script is being run from the _module directory_, which is one level up
#  from this script's directory.
# This ensures access to the module renv environment
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}
cd ..

# Define directories
notebook_template_dir="cnv-exploratory-notebooks"

# Define test data string to use with 06_infercnv.R
if [[ $TESTING -eq 1 ]]; then
  test_string="--testing"
else
  test_string=""
fi

for sample_id in SCPCS000179 SCPCS000184 SCPCS000194 SCPCS000205 SCPCS000208; do

  # define notebook output directory
  output_dir="${notebook_template_dir}/${sample_id}"
  mkdir -p ${output_dir}

  ##############################################################################
  #################################`copykat`  ##################################
  ##############################################################################


  # We run and explore copykat using euclidian distance parameter and normal cell as reference
  Rscript scripts/05_copyKAT.R --sample_id ${sample_id} --n_core ${THREADS} --distance "euclidean" --use_reference "ref"

  # We run and explore copykat using spearman distance parameter and normal cell as reference
  Rscript scripts/05_copyKAT.R --sample_id ${sample_id} --n_core ${THREADS} --distance "spearman" --use_reference "ref"

  # We run and explore copykat using euclidian distance parameter and without normal cell as reference
  Rscript scripts/05_copyKAT.R --sample_id ${sample_id} --n_core ${THREADS} --distance "euclidean" --use_reference "noref"

  # We run and explore copykat using spearman distance parameter and without normal cell as reference
  Rscript scripts/05_copyKAT.R --sample_id ${sample_id} --n_core ${THREADS} --distance "spearman" --use_reference "noref"

  # We explore `copykat` results for this sample
  Rscript -e "rmarkdown::render(input = '${notebook_template_dir}/05_copykat_exploration.Rmd',
                    params = list(sample_id = '${sample_id}', seed = 12345),
                    output_format = 'html_document',
                    output_file = '05_copykat_exploration_${sample_id}.html',
                    output_dir = '${output_dir}')"


  ##############################################################################
  #################################`infercnv`  #################################
  ##############################################################################

  # We run and explore infercnv using immune cells as reference and no HMM model
  Rscript scripts/06_infercnv.R --sample_id ${sample_id} --reference "immune" --HMM "no" ${test_string}

  # We run and explore infercnv using endothelial cells as reference and no HMM model
  Rscript scripts/06_infercnv.R --sample_id ${sample_id} --reference "endothelium" --HMM "no" ${test_string}

  # We run and explore infercnv using no normal reference and no HMM model
  Rscript scripts/06_infercnv.R --sample_id ${sample_id} --reference "none" --HMM "no" ${test_string}

  # We run and explore infercnv using both endothelial and immune cells as reference and no HMM model
  Rscript scripts/06_infercnv.R --sample_id ${sample_id} --reference "both" --HMM "no" ${test_string}

  # We run and explore infercnv using both endothelial and immune cells as reference and i3 HMM model
  Rscript scripts/06_infercnv.R --sample_id ${sample_id} --reference "both" --HMM "i3" ${test_string}

  # We run and explore infercnv using both endothelial and immune cells as reference and i6 HMM model
  Rscript scripts/06_infercnv.R --sample_id ${sample_id} --reference "both" --HMM "i6" ${test_string}

  # We run and explore infercnv using both endothelial and immune cells from all non-treated Wilms tumor patients as reference and i3 HMM model
  Rscript scripts/06_infercnv.R --sample_id ${sample_id} --reference "pull" --HMM "i3" ${test_string}


  # We explore `inferCNV` results for this sample
  Rscript -e "rmarkdown::render(input = '${notebook_template_dir}/06_cnv_infercnv_exploration.Rmd',
                    params = list(sample_id = '${sample_id}', seed = 12345),
                    output_format = 'html_document',
                    output_file = '06_infercnv_exploration_${sample_id}.html',
                    output_dir = '${output_dir}')"
done

