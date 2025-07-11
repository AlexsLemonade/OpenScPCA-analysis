#!/bin/bash
#
# Generate results and render the template exploratory notebook `reference-aggregation.Rmd`.
# This analysis only considers a subset of samples.
# They were chosen to have a range in number of cells, diversity of consensus cell type annotations, and diversity in tech/sequencing units.
# SCPCS000104: nucleus, 10Xv2
# SCPCS000105: cell, 10Xv2
# SCPCS000108: cell, 10Xv3
# SCPCS000109: cell, 10Xv3
### note that SCPCS000109 has two libraries but we only consider one of them (skip SCPCL000127)
# SCPCS000111: nucleus, 10Xv3.1
# SCPCS000115: nucleus, 10Xv3.1


set -euo pipefail

# Ensure script is being run from module root
# which is two directories up from this script's location
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}/../..

scratch_dir="scratch"
results_dir="results"
notebook_dir="exploratory-notebooks/singler-testing"
html_dir="${notebook_dir}/reference-aggregation_htmls"

# Define sample ids to process as a single string
sample_ids="SCPCS000104 SCPCS000105 SCPCS000108 SCPCS000109 SCPCS000111 SCPCS000115"

# Process with aggregated reference (default)
singler_dir="${results_dir}/singler-test_aggregated"
singler_model="${scratch_dir}/singler-model_nbatlas_aggregated.rds"
sample_ids=$sample_ids singler_results_dir="${singler_dir}" singler_model_file="${singler_model}" aggregate_singler=1 bash run-analysis.sh

# Process with non-aggregated reference
singler_dir="${results_dir}/singler-test_not-aggregated"
singler_model="${scratch_dir}/singler-model_nbatlas_not-aggregated.rds"
sample_ids=$sample_ids singler_results_dir="${singler_dir}" singler_model_file="${singler_model}" aggregate_singler=0 bash run-analysis.sh

# Render the exploratory notebook for these samples
for sample_id in $sample_ids; do
    # Loop over each sample's non-aggregated annotation just so we can grab the library_id
    for tsv_file in "${singler_dir}/${sample_id}"/*annotations.tsv; do
        library_id=$(basename "$tsv_file" | sed 's/_singler-annotations.tsv$//')

        # Skip the library from SCPCS000109 we're not considering here
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
