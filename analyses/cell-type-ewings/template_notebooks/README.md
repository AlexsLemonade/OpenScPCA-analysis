# Template notebooks

This folder contains any template notebooks that are rendered as part of a workflow in this module.

1. `cnv-workflow`: This folder contains all template notebooks used in `cnv-annotation.sh`.
2. `auc-singler-workflow`: This folder contains all template notebooks used in `auc-singler-annotation.sh`.
3. `clustering-workflow`: This folder contains all template notebooks used in `evaluate-clusters.sh`. 

## Guide notebook for assigning and evaluating cell type annotations 

The `celltype-exploration.Rmd` notebook is meant to be used as a guide for assigning and evaluating the final cell type annotations for each library in `SCPCP000015`. 
Instructions for using this guide: 

1. Ensure that you have a local copy of the results from `aucell-singler-annotation.sh`, `evaluate-clusters.sh` and `run-aucell-ews-signatures.sh` saved to `results`. 
2. Copy the contents of this notebook to a new notebook titled `<library_id>_celltype-exploration.Rmd` and save in `exploratory_analysis/final_annotation_notebooks`. 
3. Replace the `sample_id` and `library_id` with the correct IDs in the `params` list. 
4. Optionally, you may choose to update the choices for clustering based on the results from `evaluate-clusters.sh`. 
All clusters used will be calculated with the Leiden algorithm and the modularity objective function. 
To modify the nearest neighbors (default: 20) and resolution (default: 0.5) chosen use the `cluster_nn` and `cluster_res` params. 
5. Run through the notebook and update any sections of the notebook marked with `**Manual exploration**`. 
6. Render the completed notebook which will produce the rendered `html` file and a TSV with cell type annotations for that library. 
