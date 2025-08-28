#!/usr/bin/env Rscript

# This script is used to create the reference table used for assigning consensus cell types 
# the table will contain one row for each cell type combination between panglao and celldex
# where a consensus label was assigned 

# Paths ------------------------------------------------------------------------
repository_base <- rprojroot::find_root(rprojroot::is_git_root)
module_base <- rprojroot::find_root(rprojroot::is_renv_project)

# cell ontology ref file 
panglao_ref_file <- file.path(module_base, "references", "panglao-cell-type-ontologies.tsv")

# blueprint ref file 
blueprint_ref_file <- file.path(module_base, "references", "blueprint-mapped-ontologies.tsv")

# scimilarity ref file
scimilarity_file <- file.path(repository_base, "analyses", "cell-type-scimilarity", "references", "scimilarity-mapped-ontologies.tsv")

# output ref file 
consensus_ref_file <- file.path(module_base, "references", "consensus-cell-type-reference.tsv")

# Prep references --------------------------------------------------

# grab obo file
cl_ont <- ontologyIndex::get_ontology("http://purl.obolibrary.org/obo/cl/releases/2024-09-26/cl-basic.obo") 

# set up the graph to use for assigning LCA terms 
parent_terms <- cl_ont$parents
cl_graph <- igraph::make_graph(rbind(unlist(parent_terms), rep(names(parent_terms), lengths(parent_terms))))

# read in panglao file
panglao_df <- readr::read_tsv(panglao_ref_file) |>
  # rename columns to have panglao in them for easy joining later
  dplyr::select(
    panglao_ontology = "ontology_id",
    panglao_annotation = "human_readable_value",
    original_panglao_name = "panglao_cell_type" # keep original name since some map to the same ontology ID
  ) |> 
  # remove any cell types that don't have ontologies 
  tidyr::drop_na() 

# read in blueprint data 
blueprint_df <- readr::read_tsv(blueprint_ref_file)


# read in scimilarity 
scimilarity_df <- readr::read_tsv(scimilarity_file) |> 
  dplyr::select(
    scimilarity_ontology = "scimilarity_celltype_ontology",
    scimilarity_annotation = "cl_annotation",
    original_scimilarity_name = "scimilarity_celltype_annotation"
  )

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

# Create combinations ----------------------------------------------------------

# get a data frame with all combinations of panglao, blueprint, and scimilarity terms
# one row for each combination
all_ref_df <- expand.grid(
  panglao_ontology = panglao_df$panglao_ontology,
  blueprint_ontology = blueprint_df$blueprint_ontology,
  scimilarity_ontology = scimilarity_df$scimilarity_ontology
) |> 
  # remove any potentially redundant combinations based on ontology Ids 
  dplyr::distinct() |> 
  # create a column with all possible pairs between the three tools 
  # name with the cell type methods
  dplyr::mutate(cellassign_singler_pair = glue::glue("{panglao_ontology};{blueprint_ontology}"),
                cellassign_scimilarity_pair = glue::glue("{panglao_ontology};{scimilarity_ontology}"),
                singler_scimilarity_pair = glue::glue("{blueprint_ontology};{scimilarity_ontology}")) |> 
  tidyr::pivot_longer(cols = c(cellassign_singler_pair, cellassign_scimilarity_pair, singler_scimilarity_pair), names_to = "pair_id", values_to = "pairs")

# get all possible pairs from each combination and get lca for each pair that is unique 
lca_pairs_df <- all_ref_df |> 
  dplyr::select(pairs) |> 
  dplyr::distinct() |> 
  tidyr::separate(pairs, into = c("id1", "id2"), sep = ";") |> 
  dplyr::rowwise() |>
  dplyr::mutate(
    # least common shared ancestor
    lca = list(rownames(ontoProc::findCommonAncestors(id1, id2, g = cl_graph)))
  )

# Get LCA and descendants ------------------------------------------------------

# add lca and total number of descendants to data frame
# expand to have one row per unique combo + unique lca
# remember sometimes there are multiple LCA that will get returned 
lca_df <- lca_pairs_df |>
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
  dplyr::mutate(lca = dplyr::if_else(id1 == id2, id1, lca)) |> 
  # join in information for each of the lca terms including name and number of descendants
  dplyr::left_join(cl_df, by = c("lca" = "cl_ontology")) |> 
  # add in legible names for id1/id2 pairs 
  dplyr::mutate(
    id1_annotation = cl_ont$name[id1],
    id2_annotation = cl_ont$name[id2]
  )

# Filter LCAs ------------------------------------------------------------------

# vector of celltypes with <= 170 descendants that are too broad to keep as final consensus cell types
celltypes_to_exclude <- c(
  "bone cell", 
  "lining cell", 
  "blood cell", 
  "progenitor cell", 
  "supporting cell", 
  "biogenic amine secreting cell", 
  "protein secreting cell", 
  "extracellular matrix secreting cell"
  )

# create a table with all pairs and their LCA
filtered_lca_df <- lca_df |> 
  # first set criteria to keep cell types, > 170, celltypes of interest 
  dplyr::filter(total_descendants <= 170 & !(cl_annotation %in% celltypes_to_exclude) 
                | cl_annotation %in% c("neuron","columnar/cuboidal epithelial cell", "endo-epithelial cell")
                # only keep epithelial if not keratinocyte
                | (cl_annotation == "epithelial cell" & id1_annotation != "keratinocyte" & id2_annotation != "keratinocyte")) |> 
  dplyr::group_by(id1, id2) |> 
  # get the cl_annotation with the minimum value in total_descendants column 
  # discard all other rows that aren't the minimum value 
  dplyr::slice_min(total_descendants, with_ties = FALSE) |> 
  dplyr::ungroup()

# check for dups to be sure we didn't make a mistake 
dup <- filtered_lca_df |> 
  dplyr::select(id1, id2) |> 
  duplicated()

if(sum(dup) > 0) {
  stop("Duplicate LCA matches were found!")
}

# Create consensus reference ---------------------------------------------------

# create the final table with each pair all ontology and annotations and final consensus label if all three are present 

# combine lca with all combinations 
combined_ref_df <- all_ref_df |> 
  # split pairs to join with final lca assignment
  tidyr::separate(col = "pairs", into = c("id1", "id2"), sep = ";") |> 
  dplyr::left_join(filtered_lca_df, by = c("id1", "id2")) |> 
  # get rid of any combinations that don't have any possible matches 
  tidyr::drop_na()

# first get a dataframe of the minimum descendants, keeping all ties 
consensus_label_df <- combined_ref_df |> 
  # first get rid of any cases where more than one pair gives the same lca term 
  dplyr::select(panglao_ontology, blueprint_ontology, scimilarity_ontology, lca, cl_annotation, total_descendants) |> 
  dplyr::distinct() |> 
  # now select the lca term with the minimum total descendants 
  dplyr::group_by(panglao_ontology, blueprint_ontology, scimilarity_ontology) |> 
  # set a tiebreaker when HSPC is tied with stem cell
  dplyr::mutate(
    tie_priority = dplyr::if_else(cl_annotation == "hematopoietic precursor cell", 0, 1)
  ) |> 
  dplyr::slice_min(order_by = dplyr::across(c(total_descendants, tie_priority)), n = 1) |> 
  dplyr::ungroup() |> 
  dplyr::rename(
    "consensus_annotation" = cl_annotation,
    "consensus_ontology" = lca
  ) |> 
  dplyr::select(-c(tie_priority, total_descendants))

# reformat final table so that we have annotation and cl id for each method 
# we also will include the lca for each pair so if just 2 of 3 methods are present then that will be used instead of the consensus annotation 
final_consensus_df <- combined_ref_df |> 
  dplyr::select(ends_with("_ontology"), pair_id, ontology = "lca", annotation = "cl_annotation") |> 
  tidyr::pivot_wider(
    names_from = pair_id,
    values_from = c(ontology, annotation),
    # make sure the pair id is before annotation or ontology 
    names_glue = "{pair_id}_{.value}"
  ) |> 
  dplyr::left_join(consensus_label_df, by = c("panglao_ontology", "blueprint_ontology", "scimilarity_ontology"))

# final check for dups to be sure we didn't make a mistake 
dup <- final_consensus_df |> 
  dplyr::select(panglao_ontology, blueprint_ontology, scimilarity_ontology) |> 
  duplicated()

if(sum(dup)) {
  stop("Duplicate cell type combinations were found after assigning consensus cell type labels!")
}

# export table
readr::write_tsv(final_consensus_df, consensus_ref_file)
