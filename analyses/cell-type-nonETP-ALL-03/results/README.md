# Results directory instructions

Run the following codes to sync the results in a lightsail virtual computer:

```         
cd /home/lightsail-user/OpenScPCA-analysis
scripts/sync-results.py cell-type-nonETP-ALL-03 --bucket researcher-650251722463-us-east-2
```

The final rds files for this module will be found in the S3 bucket, `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/results`, and the plots can be found in `s3://researcher-650251722463-us-east-2/cell-type-nonETP-ALL-03/plots/00-01_processing_rds`.
