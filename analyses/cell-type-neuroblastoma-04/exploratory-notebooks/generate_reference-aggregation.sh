#!/bin/bash
#
# Generate results and render the template exploratory notebook `01_reference-aggregation.Rmd`.
# This analysis only considers a subset of samples chosen their range in
# number of cells, diversity of consensus cell type annotations, and diversity in tech/sequencing units.
# SCPCS000104: nucleus, 10Xv2
# SCPCS000105: cell, 10Xv2
# SCPCS000108: cell, 10Xv3
# SCPCS000109: cell, 10Xv3
### note that SCPCS000109 has two libraries but we only consider one of them (skip SCPCL000127)
# SCPCS000111: nucleus, 10Xv3.1
# SCPCS000115: nucleus, 10Xv3.1


set -euo pipefail

# Ensure script is being run from module root
# which is one directory up from this script's location
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}/..

singler_result_dir="results/singler"
notebook_dir="exploratory-notebooks"
html_dir="${notebook_dir}/reference-aggregation_htmls"

# Define sample ids to process as a single string
sample_ids="SCPCS000104 SCPCS000105 SCPCS000108 SCPCS000109 SCPCS000111 SCPCS000115"

# Process with aggregated reference (default)
sample_ids=$sample_ids bash run-analysis.sh

# Process with non-aggregated reference
sample_ids=$sample_ids aggregate_singler=0 bash run-analysis.sh

# Render the exploratory notebook for these samples
for sample_id in $sample_ids; do
    # Loop over each sample's aggregated annotation just so we can grab the library_id
    for tsv_file in "${singler_result_dir}/${sample_id}"/aggregated/*annotations.tsv; do
        library_id=$(basename "$tsv_file" | sed 's/_singler-annotations.tsv$//')

        # Skip the SCPCS000109 we're not considering here library
        if [[ $library_id == "SCPCL000127" ]]; then
            continue
        fi
        echo "Rendering exploratory notebook for sample ${sample_id} and library ${library_id}"

       Rscript -e "rmarkdown::render('${notebook_dir}/reference-aggregation.Rmd',
                                  output_file = '${library_id}_reference-aggregation.nb.html',
                                  output_dir = '${html_dir}',
                                  params = list(sample_id  = '${sample_id}',
                                                library_id = '${library_id}') )"
    done
done
