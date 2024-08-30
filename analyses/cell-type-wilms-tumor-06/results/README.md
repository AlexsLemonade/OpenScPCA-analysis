# Results directory instructions

Files in the results directory should not be directly committed to the repository.

Instead, copy results files to an S3 bucket and add a link to the S3 location in this README file.


# Azimuth compatible fetal kidney reference

To perform label transfer using Azimuth and the fetal kidney atlas, a reference is built via [`scripts/download-and-create-fetal-kidney-ref.R`](../scripts/download-and-create-fetal-kidney-ref.R) using the fetal_full.Rds object download from:
"https://datasets.cellxgene.cziscience.com/40ebb8e4-1a25-4a33-b8ff-02d1156e4e9b.rds"

The output can be found in the S3 bucket name `researcher-008971640512-us-east-2` named `references/ref.Rds` and `references/idx.annoy`.


# Clustering and label transfer from fetal references

