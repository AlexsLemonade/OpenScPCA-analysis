This directory contains files with reference information used in the module

* `nbatlas-marker-genes.tsv`: TSV file with marker genes for cell types validation taken from the `NBAtlas` paper (<https://doi.org/10.1016/j.celrep.2024.114804>).
This file was created with `../scripts/setup/prepare-nbatlas-gene-lists.R`
* `nbatlas-label-map.tsv`: TSV file mapping corresponding NBAtlas reference labels at different levels of organization between lower-level and broader family groupings, where relevant, to support cell type annotation and visualization.
  * Families considered include `cDC`, `Monocytes`, `T cells`, and `NK cells`; for cell types which do not belong to a broader family, we specify the label
  * If a cell type does not belong to a family, we use the label itself as its family
  * Broader family groupings are also mapped to themselves to support visualization code
This file was manually created based on information in the NBAtlas manuscript
* `nbatlas-ontology-ids.tsv`: TSV file with presumed Cell Ontology (`CL`) ontology IDs for labels in the `NBAtlas` reference.
The associated `CL` annotation itself is also included.
This file was manually created as described below

## Ontology ids

This section describes certain decisions made while creating `nbatlas-ontology-ids.tsv`.

* The following cell types were not assigned ontology ids but labeled `NA`:
  * `Immune cycling`: There is no equivalent in `CL`
  * `Stromal other` : The `NBAtlas` manuscript describes these cells as a small subset of stromal cells expressing adrenal cortex markers, which is a much more narrow definition than the broader `CL` ontology [`stromal cell`](http://purl.obolibrary.org/obo/CL_0000499).
* We made the following choices for these ids in particular which did not always have unambiguous `CL` id equivalents:
  * We assigned `cDC2/DC3` the id `CL:0000451`, which corresponds to `dendritic cell`.
  Because `cDC2` and `DC3` come from different developmental lineages in spite of their transcriptional similarity, we use this ontology id to represent this cell type at a very broad level. (Sources: <https://doi.org/10.1016/j.celrep.2024.114804>; <https://doi.org/10.1016/j.immuni.2023.07.001>; <https://doi.org/10.1016/j.immuni.2024.05.007>)
  * We assigned all natural killer cell types the same id `CL:0000623`, which corresponds to `natural killer cell`
  * We assigned `NKT cell` the id `CL:0000814`, which corresponds specifically to `mature NK T cell`
  * We assigned, respectively, `CD4+ T cell` and `CD8+ T cell` to `CL:0000624` and `CL:0000625` which correspond specifically to `CD4-positive, alpha-beta T cell` and `CD8-positive, alpha-beta T cell`; as such, we assumed these are `alpha-beta` (not `gamma`) CD4/8+ T cells