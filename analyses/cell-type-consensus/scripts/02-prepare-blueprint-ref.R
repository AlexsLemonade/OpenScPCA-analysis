#!/usr/bin/env Rscript

# This script is used to match cell type ontology IDs present in BlueprintEncodeData to their associated name in CL
# The output is a TSV file that contains `blueprint_ontology` and `blueprint_annotation_cl`

# paths ------------------------------------------------------------------------
module_base <- rprojroot::find_root(rprojroot::is_renv_project)

# output ref file 
blueprint_ref_file <- file.path(module_base, "references", "blueprint-mapped-ontologies.tsv")


# Map ontology ids -------------------------------------------------------------

# grab blueprint from celldex
blueprint_ref <- celldex::BlueprintEncodeData()

# grab obo file, we need this to map the ontologies from blueprint
cl_ont <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl/releases/2024-09-26/cl-basic.obo")

# get ontologies and human readable name into data frame for blueprint
# in scpca-nf we don't include the cl name so this lets us add it in
blueprint_df <- data.frame(
  blueprint_ontology = blueprint_ref$label.ont,
  blueprint_annotation_cl = cl_ont$name[blueprint_ref$label.ont]
) |>
  unique() |>
  tidyr::drop_na()

# save to output file 
readr::write_tsv(blueprint_df, blueprint_ref_file)
