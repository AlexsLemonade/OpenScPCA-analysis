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

| Field                          | Contents                        | Simulation plan                         |
| ------------------------------ | ------------------------------- | --------------------------------------- |
| `library_id`                   | string                          | keep                                    |
| `sample_id`                    | string                          | keep                                    |
| `project_id`                   | string                          | keep                                    |
| `salmon_version`               | string                          | keep                                    |
| `reference_index`              | string                          | keep                                    |
| `total_reads`                  | integer                         | keep (recalculation data not available) |
| `mapped_reads`                 | integer                         | keep (recalculation data not available) |
| `mapping_tool`                 | string                          | keep                                    |
| `alevinfry_version`            | string                          | keep                                    |
| `af_permit_type`               | string                          | keep                                    |
| `af_resolution`                | string                          | keep                                    |
| `usa_mode`                     | boolean                         | keep                                    |
| `af_num_cells`                 | string                          | keep                                    |
| `tech_version`                 | string                          | keep                                    |
| `assay_ontology_term_id`       | string                          | keep                                    |
| `seq_unit`                     | string                          | keep                                    |
| `transcript_type`              | vector (string)                 | keep                                    |
| `include_unspliced`            | boolean                         | keep                                    |
| `sample_metadata`              | data frame of metadata          | mostly keep, see below                  |
| `sample_type`                  | string                          | keep                                    |
| `filtering_method`             | string                          | keep                                    |
| `miQC_model`                   | `flexmix` model object          | remove                                  |
| `prob_compromised_cutoff`      | float                           | keep                                    |
| `scpca_filter_method`          | string                          | keep                                    |
| `min_gene_cutoff`              | integer                         | keep                                    |
| `normalization`                | string                          | keep                                    |
| `highly_variable_genes`        | vector(string)                  | keep                                    |
| `cluster_algorithm`            | string                          | keep                                    |
| `cluster_weighting`            | string                          | keep                                    |
| `cluster_nn`                   | integer                         | keep                                    |
| `singler_results`              | `DataFrame` of SingleR results  | reduce rows to simulated cell count     |
| `singler_reference`            | string                          | keep                                    |
| `singler_reference_label`      | string                          | keep                                    |
| `singler_reference_source`     | string                          | keep                                    |
| `singler_reference_version`    | string                          | keep                                    |
| `celltype_methods`             | vector(string)                  | keep                                    |
| `cellassign_predictions`       | dataframe of CellAssign results | reduce rows to simulated cell count     |
| `cellassign_reference`         | string                          | keep                                    |
| `cellassign_reference_source`  | string                          | keep                                    |
| `cellassign_reference_version` | string                          | keep                                    |
| `cellassign_reference_organs`  | string                          | keep                                    |
|                                |                                 |                                         |

### Fields in `metadata(sce)$sample_metadata`

| field                                      | contents | simulation plan |
| ------------------------------------------ | -------- | --------------- |
| `sample_id`                                | string   | keep            |
| `scpca_project_id`                         | string   | keep            |
| `submitter_id`                             | string   | keep            |
| `participant_id`                           | string   | empty string    |
| `submitter`                                | string   | keep            |
| `age`                                      | float    | keep            |
| `sex`                                      | string   | keep            |
| `diagnosis`                                | string   | keep            |
| `subdiagnosis`                             | string   | keep            |
| `tissue_location`                          | string   | keep            |
| `disease_timing`                           | string   | keep            |
| `organism`                                 | string   | keep            |
| `is_xenograft`                             | boolean  | keep            |
| `is_cell_line`                             | boolean  | keep            |
| `development_stage_ontology_term_id`       | string   | keep            |
| `sex_ontology_term_id`                     | string   | keep            |
| `organism_ontology_id`                     | string   | keep            |
| `self_reported_ethnicity_ontology_term_id` | string   | keep            |
| `disease_ontology_term_id`                 | string   | keep            |
| `tissue_ontology_term_id`                  | string   | keep            |
| `library_id`                               | string   | keep            |

## Fields in `colData(sce)`

| field                            | contents | simulation plan                 |
| -------------------------------- | -------- | ------------------------------- |
| `barcodes`                       | string   | keep                            |
| `sum`                            | float    | recalculate                     |
| `detected`                       | integer  | recalculate                     |
| `subsets_mito_sum`               | float    | calculate using percent of sum  |
| `subsets_mito_detected`          | integer  | keep                            |
| `subsets_mito_percent`           | float    | keep                            |
| `altexps_*_sum`                  | float    | recalculate                     |
| `altexps_*_detected`             | integer  | recalculate                     |
| `altexps_*_percent`              | float    | recalculate                     |
| `total`                          | float    | recalculate                     |
| `prob_compromised`               | float    | keep                            |
| `miQC_pass`                      | string   | keep                            |
| `scpca_filter`                   | string   | keep                            |
| `sizeFactor`                     | float    | keep                            |
| `submitter_celltype_annotation`  | string   | randomize from available levels |
| `cluster`                        | factor   | randomize from available levels |
| `singler_celltype_ontology`      | string   | randomize from available levels |
| `singler_celltype_annotation`    | string   | match ontology                  |
| `cellassign_celltype_annotation` | string   | randomize from available levels |
| `cellassign_max_prediction`      | float    | keep                            |

## Fields in `rowData(sce)`

| field         | contents | simulation plan |
| ------------- | -------- | --------------- |
| "gene_ids"    | string   | keep            |
| "gene_symbol" | string   | keep            |
| "mean"        | float    | recalculate     |
| "detected"    | float    | recalculate     |

## Fields in `metadata(altExp(sce))`

Similar to `metadata(sce)`, but with no sample-related data beyond SCPCX ids.
Mostly mapping-related information and counts.
`ambient_profile` contains a list of ambient ADT levels.

All of these fields are left unchanged.

## Fields in `rowData(altExp(sce))`

| field         | contents | simulation plan |
| ------------- | -------- | --------------- |
| "adt_id"      | string   | keep            |
| "target_type" | string   | keep            |
| "mean"        | float    | recalculate     |
| "detected"    | float    | recalculate     |
