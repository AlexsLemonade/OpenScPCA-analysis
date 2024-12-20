# References

This folder contains all reference files used for generating consensus cell types.

1. `panglao-cell-type-ontologies.tsv`: This file contains a table with all possible cell types in the reference used when running `CellAssign`.
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
