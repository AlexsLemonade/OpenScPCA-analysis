# Scripts 

This folder contains all scripts used for generating consensus cell types. 

1. `00-download-panglao-ref.sh`: This script is used to download the `PanglaoDB` reference file from `AlexsLemonade/scpca-nf`. 
This reference file was originally obtained from `PanglaoDB` and contains a table with all marker genes for all cell types that were used to build the references used when running `CellAssign`. 
The file will be stored in `references/PanglaoDB_markers_2020-03-27.tsv`. 

2. `01-prepare-cell-type-ontologies.R`: This script is used to assign [cell type ontologies](https://www.ebi.ac.uk/ols4/ontologies/cl) to cell types in the `PanglaoDB` reference file. 
Any cell types whose human readable label matches the value in the `cell type` column of the reference file (downloaded using the `00-download-panglao-ref.sh` file) are programmatically assigned. 
Ontology terms and labels along with the `cell type` label from the reference file are saved to a new file, `references/panglao-cell-type-ontologies.tsv`. 

3. `02-prepare-consensus-reference.R`: This script is used to create a table with all consensus cell types. 
The output table will contain one row for each combination of cell types in `PanglaoDB` and `BlueprintEncodeData` from `celldex` where a consensus cell type was identified.  
If the combination is not included in the reference file, then no consensus cell type is assigned and can be set to "Unknown". 

4. `03-save-coldata.R`: This script is used to grab the cell type annotations from the `colData` of an individual processed SCE object and save the output to a TSV file. 

5. `04-combine-celltype-tables.R`: This script is used to combine individual TSV files with cell type annotations (output by `03-save-coldata.R`) into a single TSV file. 
The consensus cell type reference is used to assign consensus cell types to all cells in the combined data frame and saved in the output TSV file. 
