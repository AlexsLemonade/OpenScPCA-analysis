# Results directory instructions

Run the following codes to sync the results in a Lightsail virtual computer:

```         
cd /home/lightsail-user/OpenScPCA-analysis
scripts/sync-results.py cell-type-nonETP-ALL-03 --bucket researcher-650251722463-us-east-2
```

These are the generated outputs for each sample in the S3 bucket:

-   `rds` objects: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results/rds`
-   metadata and ScType results: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results/`
-   umap and dot plots: `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/plots`
