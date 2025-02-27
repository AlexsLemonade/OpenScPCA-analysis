# Results directory instructions

Files in the results directory should not be directly committed to the repository.

Instead, copy results files to an S3 bucket and add a link to the S3 location in this README file.

- On a Lightsail virtual computer, run following script to sync the results:
```bash
export OPENSCPCA_RESULTS_BUCKET=researcher-009160072044-us-east-2
cd /home/lightsail-user/repo/OpenScPCA-analysis
scripts/sync-results.py cell-type-wilms-tumor-14 \
    --bucket ${OPENSCPCA_RESULTS_BUCKET}
```
#### 00. Pre-processing the provided SCE objects
No result files.

#### 01. Anchor transfer using Seurat
* The label transfer analysis was performed in two levels: `celltype` and `compartment`.
* Results are uploaded to `s3://researcher-009160072044-us-east-2/cell-type-wilms-tumor-14/results/01_anchor_transfer_seurat`


  * `[sample_id]_[level].csv` label transfer result table including cell ID, predicted cell type, along with predicted scores.
* Plots are uploaded to `s3://researcher-009160072044-us-east-2/cell-type-wilms-tumor-14/plots/01_anchor_transfer_seurat/`
  * `[sample_id]_[level].pdf` label transfer result plots consisting of 3 pages:
    1. UMAP visualization colored by transferred labels and Seurat clusters, as well as a bar plot showing cell type composition of each Seurat cluster.
    2. UMAP visualization colored and split by transferred labels.
    3. Distribution for max prediction score. Note: predictions with scores < 0.5 would be labeled as "Unknown" in this analysis.

* Anchor transfer was performed with two normalization methods in subfolders:
  * `results/01_anchor_transfer_seurat/RNA`: Results generated by normalization method `LogNormalize`.
  * `results/01_anchor_transfer_seurat/SCT`: Results generated by normalization method `SCTransform`.

* Anchor transfer results for all 10 samples are available on S3 `s3://researcher-009160072044-us-east-2/cell-type-wilms-tumor-14/results/01_anchor_transfer_seurat/summary_results.csv`.
  * Rough tumor cell classification is inferred based on anchor transfer results:
    - Cells annotated as "Cap mesenchyme" were labeled as tumor cells: CM should have gone since born [(Ref)](https://pmc.ncbi.nlm.nih.gov/articles/PMC4675508/).
    - Cells annotated as "Immune" or "Endothelium" were labeled as normal cells: Wilms tumor shows a triphasic type, meaning blastemal, stromal, and epithelial cells can be tumor cells ([Ref1](https://pmc.ncbi.nlm.nih.gov/articles/PMC6076422/) and [Ref2](https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2023.1137346/full)). In contrast, "Immune" or "Endothelium" cells are regarded as normal cells.
    - Other cells were labeled as "unknown".

