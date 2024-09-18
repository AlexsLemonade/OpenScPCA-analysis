# Results directory instructions

Files in the results directory should not be directly committed to the repository.

Instead, copy results files to an S3 bucket and add a link to the S3 location in this README file.

- On a Lightsail virtual computer, run following script to sync the results:
```bash
export OPENSCPCA_RESULTS_BUCKET=researcher-009160072044-us-east-2
cd /home/lightsail-user/git/OpenScPCA-analysis
scripts/sync-results.py cell-type-wilms-tumor-14 \
    --bucket ${OPENSCPCA_RESULTS_BUCKET}
```
#### 00. Pre-processing the provided SCE objects
No result files.

#### 01. Anchor transfer using Seurat
TBD

#### 02. Curating marker gene lists
TBD

#### 03. Cell type annotation with marker gene lists
TBD

#### 04. Tumor cell identification
TBD

#### 05. Sample merging and validation
TBD