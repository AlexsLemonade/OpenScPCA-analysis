#!/usr/bin/env Rscript

# This script is used to create the reference table used for assigning consensus cell types 
# the table will contain one row for each cell type combination between panglao and celldex
# where a consensus label was assigned 

# Paths ------------------------------------------------------------------------
module_base <- rprojroot::find_root(rprojroot::is_renv_project)

# cell ontology ref file 
panglao_ref_file <- file.path(module_base, "references", "panglao-cell-type-ontologies.tsv")

# output ref file 
consensus_ref_file <- file.path(module_base, "references", "consensus-cell-type-reference.tsv")

# Prep references --------------------------------------------------

# grab obo file
cl_ont <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl-basic.obo") 

# set up the graph to use for assigning LCA terms 
parent_terms <- cl_ont$parents
cl_graph <- igraph::make_graph(rbind(unlist(parent_terms), rep(names(parent_terms), lengths(parent_terms))))

# read in panglao file
panglao_df <- readr::read_tsv(panglao_ref_file) |>
  # rename columns to have panglao in them for easy joining later
  dplyr::select(
    panglao_ontology = "ontology_id",
    panglao_annotation = "human_readable_value"
  ) |> 
  # remove any cell types that don't have ontologies 
  tidyr::drop_na()

# grab singler ref from celldex
blueprint_ref <- celldex::BlueprintEncodeData()

# get ontologies and human readable name into data frame
blueprint_df <- data.frame(
  blueprint_ontology = blueprint_ref$label.ont,
  blueprint_annotation_main = blueprint_ref$label.main,
  blueprint_annotation_fine = blueprint_ref$label.fine
) |>
  unique() |> 
  tidyr::drop_na()

# Get LCA and descendants ------------------------------------------------------

# get total descendants for each term in CL  
# turn cl_ont into data frame with one row per term
cl_df <- data.frame(
  cl_ontology = cl_ont$id,
  cl_annotation = cl_ont$name
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    descendants = list(ontologyIndex::get_descendants(cl_ont, cl_ontology, exclude_roots = TRUE)),
    total_descendants = length(descendants)
  )

# get a data frame with all combinations of panglao and blueprint terms
# one row for each combination
all_ref_df <- expand.grid(
  panglao_df$panglao_ontology,
  blueprint_df$blueprint_ontology
) |>
  dplyr::rename(
    panglao_ontology = "Var1",
    blueprint_ontology = "Var2"
  ) |>
  # add in the human readable values for each ontology term
  # account for ontologies showing up multiple times 
  dplyr::left_join(blueprint_df, by = "blueprint_ontology", relationship = "many-to-many") |>
  dplyr::left_join(panglao_df, by = "panglao_ontology", relationship = "many-to-many") |> 
  unique() # only keep unique combinations

# add lca and total number of descendants to data frame
# expand to have one row per unique combo + unique lca
# later we will remove the extra lca assignments so that there is only one row per combination and one consensus label 
lca_df <- all_ref_df |>
  dplyr::rowwise() |>
  dplyr::mutate(
    # least common shared ancestor
    lca = list(rownames(ontoProc::findCommonAncestors(blueprint_ontology, panglao_ontology, g = cl_graph)))
  ) |> 
  dplyr::mutate(
    total_lca = length(lca), # get total number for filtering later 
    lca = paste0(lca, collapse = ",") # make it easier to split the lca terms
  ) |>
  # split each lca term into its own column
  tidyr::separate(lca, into = c("lca_1", "lca_2", "lca_3"), sep = ",") |>
  # transpose so that instead of lca being in a column there is one row per lca 
  tidyr::pivot_longer(
    cols = dplyr::starts_with("lca"),
    names_to = "lca_number",
    values_to = "lca"
  ) |>
  tidyr::drop_na() |>
  dplyr::select(-lca_number) |> # don't need this column 
  # account for any cases where the ontology IDs are exact matches
  # r complains about doing this earlier since the lca column holds lists until now
  dplyr::mutate(lca = dplyr::if_else(blueprint_ontology == panglao_ontology, blueprint_ontology, lca)) |> 
  # join in information for each of the lca terms including name and number of descendants
  dplyr::left_join(cl_df, by = c("lca" = "cl_ontology"))

# Set consensus labels ---------------------------------------------------------

# get a table with only combinations that will have an assigned consensus label 
consensus_labels_df <- lca_df |> 
  # everything with more than 1 lca gets removed with the exception of HSCs
  dplyr::filter(total_lca <=1 | cl_annotation == "hematopoietic precursor cell") |> 
  # keep everything with total descendants < 170 except for neuron and epithelial cell when blueprint calls it as epithelial 
  dplyr::filter(total_descendants <= 170 | cl_annotation %in% c("neuron", "epithelial cell") & blueprint_annotation_main == "Epithelial cells") |> 
  # get rid of terms that have low number of descendants but are still too broad 
  dplyr::filter(!(cl_annotation %in% c("bone cell", "lining cell", "blood cell", "progenitor cell", "supporting cell"))) |> 
  dplyr::select(panglao_ontology, blueprint_ontology, consensus_ontology = lca, consensus_annotation = cl_annotation) 

# export table
readr::write_tsv(consensus_labels_df, consensus_ref_file)
