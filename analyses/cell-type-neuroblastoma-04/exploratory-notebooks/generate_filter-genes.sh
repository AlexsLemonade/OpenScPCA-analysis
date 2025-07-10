#!/bin/bash
#
# Generate results and render the template exploratory notebook `filter-genes.Rmd`.
# This analysis only considers a subset of samples.
# They were chosen to have a range in number of cells, diversity of consensus cell type annotations, and diversity in tech/sequencing units.
# SCPCS000104: nucleus, 10Xv2
# SCPCS000105: cell, 10Xv2
# SCPCS000108: cell, 10Xv3
# SCPCS000109: cell, 10Xv3
### note that SCPCS000109 has two libraries
# SCPCS000111: nucleus, 10Xv3.1
# SCPCS000115: nucleus, 10Xv3.1


set -euo pipefail

# Ensure script is being run from module root
# which is one directory up from this script's location
module_dir=$(dirname "${BASH_SOURCE[0]}")
cd ${module_dir}/..

results_dir="results"
scratch_dir="scratch"
notebook_dir="exploratory-notebooks"
html_dir="${notebook_dir}/filter-genes_htmls"
mkdir -p ${html_dir}

# Define sample ids to process as a single string
sample_ids="SCPCS000104 SCPCS000105 SCPCS000108 SCPCS000109 SCPCS000111 SCPCS000115"

# Process with reference where ribosomal and mitochondrial genes have been filtered out before SingleR training
singler_dir="${results_dir}/singler-test_filter-genes-TRUE"
singler_model="${scratch_dir}/singler-model_nbatlas_filter-genes-TRUE.rds"
sample_ids=$sample_ids singler_results_dir="${singler_dir}" singler_model_file="${singler_model}" filter_genes_singler=1 bash run-analysis.sh

# Process with reference where no genes are filtered out before SingleR training (default)
singler_dir="${results_dir}/singler-test_filter-genes-FALSE"
singler_model="${scratch_dir}/singler-model_nbatlas_filter-genes-FALSE.rds"
sample_ids=$sample_ids singler_results_dir="${singler_dir}" singler_model_file="${singler_model}" filter_genes_singler=0 bash run-analysis.sh


# Render the exploratory notebook for these samples
for sample_id in $sample_ids; do
    # Loop over each sample's annotation just so we can grab the library_id
    for tsv_file in "${singler_dir}/${sample_id}"/*annotations.tsv; do
        library_id=$(basename "$tsv_file" | sed 's/_singler-annotations.tsv$//')

        echo "Rendering exploratory notebook for sample ${sample_id} and library ${library_id}"

        Rscript -e "rmarkdown::render('${notebook_dir}/filter-genes.Rmd',
                                    output_file = '${library_id}_filter-genes.nb.html',
                                    output_dir = '${html_dir}',
                                    params = list(sample_id  = '${sample_id}',
                                                    library_id = '${library_id}') )"
    done
done
