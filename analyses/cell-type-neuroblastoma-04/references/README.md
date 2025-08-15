This directory contains files with reference information used in the module

* `nbatlas-marker-genes.tsv`: Marker genes for cell types validation taken from the `NBAtlas` paper (<https://doi.org/10.1016/j.celrep.2024.114804>).
This file was created with `../scripts/setup/prepare-nbatlas-gene-lists.R`
* `nbatlas-label-map.tsv`: TSV mapping corresponding NBAtlas reference labels at different labels of organization between the main labels and "immune zoom" labels.
This file was manually created based on information in the NBAtlas manuscript
* `nbatlas-ontology-ids.tsv`: TSV with presumed ontology IDs for labels in the `NBAtlas` reference.
This file was manually created; `NA` values represent cell types for which we could not confidently determine an exact ontology ID.
  * Note that the `Stromal other` cell type is given an `NA` ontology even though there exists an ontology ID for this cell type, because as defined in the manuscript this grouping does not represent all stromal cells but a subset of stromal cells expressing adrenal cortex markers
