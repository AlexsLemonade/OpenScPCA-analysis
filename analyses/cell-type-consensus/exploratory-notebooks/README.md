# Exploratory notebooks

This folder contains exploratory notebooks for this module. 

1. `01-reference-exploration.Rmd`: This notebook was used to explore possible consensus label assignments between cell types in the `PanglaoDB` and `BlueprintEncodeData` references. 
Observations made in this notebook were used to define the set of possible consensus labels to be included in [`references/consensus-cell-type-reference.tsv`](../references/consensus-cell-type-reference.tsv). 

2. `02-explore-consensus-results.Rmd`: This notebook summarizes the consensus labels assigned to all ScPCA samples. 
Prior to rendering this notebook results from the `cell-type-consensus` module in `OpenScPCA-nf` using the `2024-11-25` were downloaded. 

3. `03-osteosarcoma-consensus-celltypes.Rmd`: This notebook summarizes the consensus labels assigned to projects with Osteosarcoma samples. 

4. `04-brain-and-cns-consensus-celltypes.Rmd`: THis notebook summarizes the consensus labels assigned to all samples belonging to projects with Brain and CNS tumors, excluding any samples that are part of multiplexed libraries. 

## Utils 

The `utils` folder contains scripts with any functions that are used by multiple notebooks. 
