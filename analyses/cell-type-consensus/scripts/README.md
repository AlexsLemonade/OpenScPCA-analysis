# Scripts 

This folder contains all scripts used for generating consensus cell types. 

1. `00-download-panglao-ref.sh`: This script is used to download the `PanglaoDB` reference file from `AlexsLemonade/scpca-nf`. 
This reference file was originally obtained from `PanglaoDB` and contains a table with all marker genes for all cell types that were used to build the references used when running `CellAssign`. 
The file will be stored in `references/PanglaoDB_markers_2020-03-27.tsv`. 

2. `01-prepare-cell-type-ontologies.sh`: This script is used to assign [cell type ontologies](https://www.ebi.ac.uk/ols4/ontologies/cl) to cell types in the `PanglaoDB` reference file. 
Any cell types whose human readable label matches the value in the `cell type` column of the reference file (downloaded using the `00-download-panglao-ref.sh` file) are programmatically assigned. 
Ontology terms and labels along with the `cell type` label from the reference file are saved to a new file, `references/panglao-cell-type-ontologies.tsv`. 
