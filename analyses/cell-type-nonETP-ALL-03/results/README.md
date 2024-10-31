# Results directory instructions

Run the following codes to sync the results in a Lightsail virtual computer:

```         
cd /home/lightsail-user/OpenScPCA-analysis
scripts/sync-results.py cell-type-nonETP-ALL-03 --bucket researcher-650251722463-us-east-2
```

These are the generated outputs for each sample in the S3 bucket:

-   `rds` objects: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results/rds`
-   metadata and ScType results: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results/`
-   CopyKat results: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results/copykat_output`
-   InferCNV results: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results/infercnv_output`
-   evaluating cluster separation, stability, and purity: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results/evalClus`
-   umap and dot plots: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/plots`
-   violin and stacked bar plots for exploring the results of CopyKat prediction: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/plots/copykat_exploration`
-   ridge plots showing the ScType score for each cell type in annotated B cells from ScType, SingleR, and CellAssign, as well as the scatter plots showing the relationship between B cell ScType score and cluster purity of these cells: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/plots/sctype_exploration`
-   ridge plots showing the ScType score for each cell type in fine-tuned B cells and feature plots showing the distribution of B cell ScType score: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/plots/sctype_exploration`
-   final submission `tsv` files and `png` for cell type and/or tumor cell classification: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results/submission_table`

\*\*All the plots are also found in the repository plots/.
