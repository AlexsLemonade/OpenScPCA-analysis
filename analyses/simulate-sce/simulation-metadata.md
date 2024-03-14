# Simulation and metadata in SCE objects

This document describes what components of the SCE object are preserved and/or modified during simulations with the `simulate-sce.R` script.
Most contents are preserved in original form, but some are recalculated from the simulated data or randomized with the intent of preserving characteristics of the input data, such as the number and content of cell type annotations.

## Data matrices

| Field                  | Contents | Simulation plan                      |
| ---------------------- | -------- | ------------------------------------ |
| counts assay           |          | simulate with `splatter` (100 cells) |
| spliced assay          |          | fraction of counts                   |
| PCA matrix             |          | recalculate (10 dimensions)          |
| UMAP matrix            |          | recalculate                          |
| altExp counts assay(s) |          | simulate with `splatter`             |

## Fields in `metadata(sce)`

| Field                        | Contents                        | Simulation plan                         |
| ---------------------------- | ------------------------------- | --------------------------------------- |
| library_id                   | string                          | keep                                    |
| sample_id                    | string                          | keep                                    |
| project_id                   | string                          | keep                                    |
| salmon_version               | string                          | keep                                    |
| reference_index              | string                          | keep                                    |
| total_reads                  | integer                         | keep (recalculation data not available) |
| mapped_reads                 | integer                         | keep (recalculation data not available) |
| mapping_tool                 | string                          | keep                                    |
| alevinfry_version            | string                          | keep                                    |
| af_permit_type               | string                          | keep                                    |
| af_resolution                | string                          | keep                                    |
| usa_mode                     | bool                            | keep                                    |
| af_num_cells                 | string                          | keep                                    |
| tech_version                 | string                          | keep                                    |
| assay_ontology_term_id       | string                          | keep                                    |
| seq_unit                     | string                          | keep                                    |
| transcript_type              | vector (string)                 | keep                                    |
| include_unspliced            | bool                            | keep                                    |
| sample_metadata              | data frame of metadata          | mostly keep, see below                  |
| sample_type                  | string                          | keep                                    |
| filtering_method             | string                          | keep                                    |
| miQC_model                   | flexmix model object            | remove                                  |
| prob_compromised_cutoff      | float                           | keep                                    |
| scpca_filter_method          | string                          | keep                                    |
| min_gene_cutoff              | integer                         | keep                                    |
| normalization                | string                          | keep                                    |
| highly_variable_genes        | vector(string)                  | keep                                    |
| cluster_algorithm            | string                          | keep                                    |
| cluster_weighting            | string                          | keep                                    |
| cluster_nn                   | integer                         | keep                                    |
| singler_results              | DataFrame of singleR statistics | reduce rows to simulated cell count     |
| singler_reference            | string                          | keep                                    |
| singler_reference_label      | string                          | keep                                    |
| singler_reference_source     | string                          | keep                                    |
| singler_reference_version    | string                          | keep                                    |
| celltype_methods             | vector(string)                  | keep                                    |
| cellassign_predictions       | dataframe of cellassign probs   | reduce rows to simulated cell count     |
| cellassign_reference         | string                          | keep                                    |
| cellassign_reference_source  | string                          | keep                                    |
| cellassign_reference_version | string                          | keep                                    |
| cellassign_reference_organs  | string                          | keep                                    |
|                              |                                 |                                         |

### Fields in `metadata(sce)$sample_metadata`

| field                                    | contents | simulation plan |
| ---------------------------------------- | -------- | --------------- |
| sample_id                                | string   | keep            |
| scpca_project_id                         | string   | keep            |
| submitter_id                             | string   | keep            |
| participant_id                           | string   | empty string    |
| submitter                                | string   | keep            |
| age                                      | float    | keep            |
| sex                                      | string   | keep            |
| diagnosis                                | string   | keep            |
| subdiagnosis                             | string   | keep            |
| tissue_location                          | string   | keep            |
| disease_timing                           | string   | keep            |
| organism                                 | string   | keep            |
| is_xenograft                             | bool     | keep            |
| is_cell_line                             | bool     | keep            |
| development_stage_ontology_term_id       | string   | keep            |
| sex_ontology_term_id                     | string   | keep            |
| organism_ontology_id                     | string   | keep            |
| self_reported_ethnicity_ontology_term_id | string   | keep            |
| disease_ontology_term_id                 | string   | keep            |
| tissue_ontology_term_id                  | string   | keep            |
| library_id                               | string   | keep            |

## Fields in `colData(sce)`

| field                          | contents | simulation plan                 |
| ------------------------------ | -------- | ------------------------------- |
| barcodes                       | string   | keep                            |
| sum                            | float    | recalculate                     |
| detected                       | integer  | recalculate                     |
| subsets_mito_sum               | float    | calculate using percent of sum  |
| subsets_mito_detected          | integer  | keep                            |
| subsets_mito_percent           | float    | keep                            |
| altexps\_\*\_sum               | float    | recalculate                     |
| altexps\_\*\_detected          | integer  | recalculate                     |
| altexps\_\*\_percent           | float    | recalculate                     |
| total                          | float    | recalculate                     |
| prob_compromised               | float    | keep                            |
| miQC_pass                      | string   | keep                            |
| scpca_filter                   | string   | keep                            |
| sizeFactor                     | float    | keep                            |
| submitter_celltype_annotation  | string   | randomize from available levels |
| cluster                        | factor   | randomize from available levels |
| singler_celltype_ontology      | string   | randomize from available levels |
| singler_celltype_annotation    | string   | match ontology                  |
| cellassign_celltype_annotation | string   | randomize from available levels |
| cellassign_max_prediction      | float    | keep                            |
