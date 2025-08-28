# Scripts 

This folder contains all scripts used for generating consensus cell types. 

1. `00-download-panglao-ref.sh`: This script is used to download the `PanglaoDB` reference file from `AlexsLemonade/scpca-nf`. 
This reference file was originally obtained from `PanglaoDB` and contains a table with all marker genes for all cell types that were used to build the references used when running `CellAssign`. 
The file will be stored in `references/PanglaoDB_markers_2020-03-27.tsv`. 

2. `01-prepare-cell-type-ontologies.R`: This script is used to assign [cell type ontologies](https://www.ebi.ac.uk/ols4/ontologies/cl) to cell types in the `PanglaoDB` reference file. 
Any cell types whose human readable label matches the value in the `cell type` column of the reference file (downloaded using the `00-download-panglao-ref.sh` file) are programmatically assigned. 
Ontology terms and labels along with the `cell type` label from the reference file are saved to a new file, `references/panglao-cell-type-ontologies.tsv`. 

3. `02-prepare-blueprint-ref.R`: This script grabs the cell types that are part of `BlueprintEncodeData` in `celldex` and saves the ontology terms and associated names from the cell type ontology (CL). 
The terms and names are saved to a new file, `references/blueprint-mapped-ontologies.tsv`.

4. `03-prepare-consensus-reference.R`: This script is used to create a table with all consensus cell types. 
The output table will contain one row for each combination of cell types in `PanglaoDB`, `BlueprintEncodeData` from `celldex`, and `SCimilarity` where a consensus cell type was identified.  
If the combination is not included in the reference file, then no consensus cell type is assigned and can be set to "Unknown". 

5. `04-assign-consensus-celltypes.R`: This script is used to grab the existing cell type annotations from the `colData` of an individual processed SCE object and assign the appropriate consensus cell type based on the `singler_celltype_ontology` (`BlueprintEncodeData`) and the `cellassign_celltype_ontology` (`PanglaoDB`). 
All annotations, including the consensus annotation, are then saved to a TSV file. 
An additional TSV file containing the gene expression for all marker genes found in `references/validation-markers.tsv` will also be saved. 

6. `00-download-cellmarker-ref.sh`: This script is used to download the [marker genes for all Human tissues from `CellMarker2.0`](http://117.50.127.228/CellMarker/CellMarker_download.html). 
This file will be stored in `references/Cell_markers_Human.xlsx`. 

7. `05-generate-validation-markers.R`: This script is used to create a table of top marker genes for each cell types represented in the consensus cell type labels. 
Prior to running this script, the marker gene file from `CellMarker2.0` must be downloaded with `00-download-cellmarker-ref.sh`. 
The output includes the top observed marker genes for all cell types in the `validation_group_ontology` column of `references/consensus-validation-groups.tsv`. 
Marker genes are ranked based on how frequently they are observed in all tissues present in `CellMarker2.0`. 
