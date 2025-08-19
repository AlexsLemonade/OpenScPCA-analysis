## Pooled results

Pooled analysis results are organized by project as follows:

```console
└── {project id}
    └── {sample id}
        ├── {first normal reference}
        │   ├── {library_id}_cnv-metadata.tsv
        │   ├── {library_id}_cnv-obj.rds
        │   ├── {library_id}_infercnv.png
        │   └── {library_id}_infercnv-results.nb.html
        └── {second normal reference}
            └── ...
```

References are named as `{cell type group}_{pooled/internal}`, where `pooled` references contain cells of the given cell type group pooled across all samples in the project, and `internal` references contain cells of that group only from the given library itself.

### Description of results files

The following files are `inferCNV` output produced by `../scripts/01_run-infercnv.R`:

* `{library_id}_cnv-metadata.tsv`
* `{library_id}_cnv-obj.rds`
* `{library_id}_infercnv.png`

The notebook files `{library_id}_infercnv-results.nb.html` were created from the given project's relevant template notebook to explore `inferCNV` results in `../template-notebooks/`.
