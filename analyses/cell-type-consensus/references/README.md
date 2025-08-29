# References

This folder contains all reference files used for generating consensus cell types.

#### PanglaoDB cell type ontology identifiers 

The `panglao-cell-type-ontologies.tsv` file contains a table with all possible cell types in the reference used when running `CellAssign`.
The table includes the following columns:

|  |   |
| --- | --- |
| `ontology_id` | [cell type (CL) ontology identifier term](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `human_readable_value` | Label associated with the cell type ontology term |
| `panglao_cell_type` | Original name for the cell type as set by `PanglaoDB` |

To generate this file follow these steps:

- Download the original reference file with `00-download-panglao-ref.sh`.
- Programmatically assign ontology labels with `01-prepare-cell-type-ontologies.sh`.
- Assign any ontology values manually by finding the most representative [cell type ontology (CL) identifier term](https://www.ebi.ac.uk/ols4/ontologies/cl).

There were a few terms that were assigned manually that did not have an obvious or exact match in the cell type ontology that should be noted:

- `Kidney progenitor cells` were assigned the [`CL:0000324` for metanephric mesenchyme stem cell](https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCL_0000324).
This term refers specifically to the stem cells that ultimately comprise the nephron, but does not account for the part of the kidney that is derived from ureteric bud cells.
- `Meningeal cells` were assigned to [`CL:0000708` for leptomeningeal cell](https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCL_0000708).
This was the closest term that covered general cell types found in the meninges.
- `Pancreatic progenitor cells` were assigned to [`CL:0002351` for progenitor cell of endocrine pancreas](https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCL_0002351).
This term only covers progenitor cells for the endocrine pancreas and does not cover the exocrine pancreas.
There were no terms that encompassed both other than `progenitor cell`.
- `Osteoclast precursor cells` were assigned to [`CL:0000576` for monocyte](https://www.ebi.ac.uk/ols4/ontologies/cl/classes/http%253A%252F%252Fpurl.obolibrary.org%252Fobo%252FCL_0000576).
Monocytes differentiate into mononuclear osteoclasts which are then activated and become multinucleated osteoclasts.
Because monocytes are the "precursor" to the differentiated osteoclast, we chose to use this term.
- `NA` was used for `Undefined placental cells` and `Transient cells` as no clear cell type from the cell ontology was identified.

#### Blueprint cell type ontology identifiers

The `blueprint-mapped-ontologies.tsv` file contains a table with all cell types from the [`BlueprintEncodeData` reference available through `celldex`](https://rdrr.io/github/LTLA/celldex/man/BlueprintEncodeData.html). 

The table includes the following columns: 

|  |   |
| --- | --- |
| `blueprint_ontology` | Cell type ontology term for `BlueprintEncodeData` cell type |
| `blueprint_annotation_cl` | Human readable value associated with the cell type ontology term for `blueprint_ontology` |

#### Consensus cell type references

The `consensus-cell-type-reference.tsv` file contains a table with all cell type combinations between the `PanglaoDB` reference (used with `CellAssign`), the `BlueprintEncodeData` reference (used with `SingleR`), and the `SCimilarity` reference (used with `SCimilarity`) for which a consensus cell type is identified. 

The table includes the following columns: 

|  |   |
| --- | --- |
| `panglao_ontology` | Cell type ontology term for `PanglaoDB` cell type |
| `panglao_annotation` | Name for the `panglao_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `original_panglao_name` | Original name for the cell type as set by `PanglaoDB` |
| `blueprint_ontology` | Cell type ontology term for `BlueprintEncodeData` cell type |
| `blueprint_annotation_cl` | Name for the `blueprint_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `scimilarity_ontology` | Cell type ontology term for `SCimilarity` cell type (obtained from [`cell-type-scimilarity/references/scimilarity-mapped-ontologies.tsv`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/analyses/cell-type-scimilarity/references/scimilarity-mapped-ontologies.tsv)) |
| `scimilarity_annotation` | Name for the `scimilarity_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |
| `original_scimilarity_name` | Original name for the cell type as set by `SCimilarity` |
| `cellassign_singler_pair_ontology` | Cell type ontology term for the consensus cell type between `PanglaoDB` (`CellAssign`) and `BlueprintEncodeData` (`SingleR`) |
| `cellassign_singler_pair_annotation` | Name for the consensus cell type between `PanglaoDB` (`CellAssign`) and `BlueprintEncodeData` (`SingleR`) |
| `cellassign_scimilarity_pair_ontology` | Cell type ontology term for the consensus cell type between `PanglaoDB` (`CellAssign`) and `SCimilarity` |
| `cellassign_scimilarity_pair_annotation` | Name for the consensus cell type between `PanglaoDB` (`CellAssign`) and `SCimilarity` |
| `singler_scimilarity_pair_ontology` | Cell type ontology term for the consensus cell type between `SCimilarity` and `BlueprintEncodeData` (`SingleR`) |
| `singler_scimilarity_pair_annotation` | Name for the consensus cell type between `SCimilarity` and `BlueprintEncodeData` (`SingleR`) |
| `consensus_ontology` | Cell type ontology term for consensus cell type |
| `consensus_annotation` | Name for the `consensus_ontology` term as defined by [Cell Ontology](https://www.ebi.ac.uk/ols4/ontologies/cl) |

This file was generated by running [`scripts/03-prepare-consensus-reference.R`](../scripts/03-prepare-consensus-reference.R). 

The `consensus-immune-cell-types.tsv` file contains all consensus cell type labels that are part of the immune component. 

The `consensus-validation-groups.tsv` file contains all possible consensus cell type labels and assigns to a broader group to use for validating cell type assignments. 

The table includes the following columns: 

|  |   |
| --- | --- |
| `consensus_ontology` | Cell type ontology term for consensus cell type |
| `consensus_annotation` | Human readable name for the consensus cell type |
| `validation_group_ontology` | Cell type ontology term for broader cell type group used for validation |
| `consensus_annotation` | Human readable name for broader cell type group used for validation |

#### Cell type validation marker genes

The `validation-markers.tsv` file contains a list of potential marker genes to aid in validating consensus cell types. 
These markers were obtained by grabbing the most frequently observed markers for each cell type included in `consensus-validation-groups.tsv` from the [`CellMarker2.0` list of human marker genes](http://117.50.127.228/CellMarker/CellMarker_download.html). 
The top 10 genes (sometimes more if there is a tie in the frequency) are included for each cell type. 

The table includes the following columns: 

|  |   |
| --- | --- |
| `validation_group_ontology` | Cell type ontology term for broader cell type group used for validation |
| `consensus_annotation` | Human readable name for broader cell type group used for validation |
| `ensembl_gene_id` | Ensembl gene identifier for the marker gene |
| `gene_symbol` | Gene symbol for the marker gene |
| `number_of_tissues` | Total number of tissues that express that marker gene in the specified cell type | 
| `celltype_total_tissues` | Total number of tissues that contained the specified cell type | 
| `percent_tissues` | Percentage of tissues that express the marker gene in the specified cell type | 

To generate this file follow these steps:

- Download the original reference file with `00-download-cellmarker-ref.sh`.
- Generate a table of marker genes with `05-generate-validation-markers.sh`.
