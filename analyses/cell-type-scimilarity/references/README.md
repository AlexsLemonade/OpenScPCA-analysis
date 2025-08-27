# References

This folder contains any reference files used for cell type annotation. 

* `scimilarity-missing-ontologies.tsv`: This file contains a table of all annotations that are part of the [`model_v1.1`](https://zenodo.org/records/10685499) that do not have a direct cell ontology match. 
For each cell type annotation present in the `scimilarity_celltype_annotation` column we identified the most closely matching name in cell ontology, indicated in the `cl_annotation` column. 
If there are new annotations found in the `SCimilarity` model that do not match cell ontology directly, they should be added to this table. 

* `scimilarity-mapped-ontologies.tsv`: This file contains a table of all annotations that are part of [`model_v1.1`](https://zenodo.org/records/10685499) and their associated [cell ontology identifiers](https://www.ebi.ac.uk/ols4/ontologies/cl).
This file was produced by running `scripts/01-assign-ontology-ids.tsv` and contains the following columns: 

| | | 
| -| - | 
| `scimilarity_celltype_annotation` | The cell type annotation present in the `SCimilarity` model | 
| `cl_annotation` | The human readable name of the cell type found in the Cell Ontology | 
| `scimilarity_celltype_ontology` | The ontology identifier associated with the name of the cell type in the Cell Ontology | 
