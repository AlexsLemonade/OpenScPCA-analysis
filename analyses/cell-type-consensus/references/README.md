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
- Programmatically assign ontology lables with `01-prepare-cell-type-ontologies.sh`. 
- Assign any ontology values manually by finding the most representive [cell type ontology (CL) identifier term](https://www.ebi.ac.uk/ols4/ontologies/cl). 
