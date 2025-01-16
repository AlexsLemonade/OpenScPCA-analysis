# Frequently asked questions

### Why didn't the sample/project I specified when running the data download script download?

First, we recommend using the `--dryrun` flag when running the [data download script](../getting-started/accessing-resources/getting-access-to-data.md#using-the-download-data-script) to check which files _would_ be downloaded.
This will confirm that there is nothing wrong with your internet connection and that you are properly [logged into your AWS profile](../technical-setup/environment-setup/configure-aws-cli.md#logging-in-to-a-new-session).

If running the script with `--dryrun` states that _only_ the `DATA_USAGE.md` file is being downloaded, this means the data files you are attempting to download do not exist.
There are two main reasons why this might occur:

1. Not all samples are available in `AnnData` format.
   If you are attempting to download `AnnData` format, please [consult the ScPCA documentation](https://scpca.readthedocs.io/en/stable/faq.html#which-samples-can-i-download-as-anndata-objects) to learn more about which types of projects do not have `AnnData` files.
1. The sample or project ID(s) that you specified during download do(es) not exist.

Therefore, you may wish to list some of the files in the data release S3 bucket to confirm the sample/project files you are trying to download actually exist.

Data files in each release are organized on S3 as:

```{ .console .no-copy title="Release file structure"}
{Release}
    ├── {Project ID}
    │   └── {Sample ID}
    │       └── {Library files}
    ├── bulk_metadata.tsv (if applicable)
    ├── bulk_quant.tsv (if applicable)
    └── single_cell_metadata.tsv
```

Follow these instructions to determine the contents of a release:

1.  First, find the names of all the data releases, which are named based on their release date in [ISO 8601 format](https://en.wikipedia.org/wiki/ISO_8601).
    You can use the [`download-data.py`](../getting-started/accessing-resources/getting-access-to-data.md#downloaded-data-file-structure) script with the `--list-releases` flag to do this.
    Don't forget to [log into your AWS `openscpca` profile first](../technical-setup/environment-setup/configure-aws-cli.md#logging-in-to-a-new-session).
    We also strongly suggest [storing your AWS profile name](../technical-setup/environment-setup/configure-aws-cli.md#storing-your-aws-profile-name) to help handle credentials.

    ```bash
    # Ensure you are in the top-level of the repository
    # cd path/to/OpenScPCA-analysis

    # If your profile name is not stored, you can set for use in this terminal session only:
    export AWS_PROFILE=openscpca

    # List all releases in the data release bucket
    ./download-data.py --list-releases
    ```

1.  Identify the most recent date from the output.
    In the example output below, the most recent release is `2024-05-01`.

    ```{ .console .no-copy title="Output from listing all releases"}
    Available release dates:
    2024-03-11
    2024-04-06
    2024-05-01
    ```

1.  All data is stored on S3 in `s3://openscpca-data-release/{release name}/`.
    To list all available files, you'll need to use the AWS CLI [in a terminal](../getting-started/project-tools/using-the-terminal.md) and run `aws s3 ls` to list the contents of a given S3 bucket.
    For example, you can list files in a given release as:

    ```bash
    # List all releases in the 2024-05-01 release
    aws s3 ls s3://openscpca-data-release/2024-05-01/
    ```

    This has the (abbreviated) output:

    ```{ .console .no-copy title="Output from listing all projects in the 2024-05-01 release"}
    PRE SCPCP000001/
    PRE SCPCP000002/
    PRE SCPCP000003/
    PRE SCPCP000004/
    ```

    Note that the `PRE` prefix appears for directories that are inside the overall data release bucket, `s3://openscpca-data-release`.

    !!! tip "Don't forget the trailing slash!"
        The trailing slash at the end of `s3://openscpca-data-release/2024-05-01/` is necessary for contents to be listed.
        If you omit the slash and just run `aws s3 ls s3://openscpca-data-release/2024-05-01`, the output will only list the directory itself, and not its contents.

1.  From there, you can continue to list the next level of nested files and directories (prefixes).
    For example, you can list all files in the `SCPCS000001` sample from the `SCPCP000001` project with:

    ```bash
    # List all SCPCS000001 files, again specifying the AWS profile
    aws s3 ls s3://openscpca-data-release/2024-05-01/SCPCP000001/SCPCS000001/
    ```

    This has the output:

    ```{ .console .no-copy title="Output from listing all SCPCS000001 files in the 2024-05-01 release"}
    2024-04-30 09:38:01    2583665 SCPCL000001_celltype-report.html
    2024-04-30 09:38:01   47346465 SCPCL000001_filtered.rds
    2024-04-30 09:40:45   48974004 SCPCL000001_filtered_rna.h5ad
    2024-04-30 09:38:01   45388563 SCPCL000001_processed.rds
    2024-04-30 09:40:45  214177811 SCPCL000001_processed_rna.h5ad
    2024-04-30 09:38:01    3057684 SCPCL000001_qc.html
    2024-04-30 09:38:03   47019247 SCPCL000001_unfiltered.rds
    2024-04-30 09:40:45   94000784 SCPCL000001_unfiltered_rna.h5ad
    ```

1.  If you attempt to list files with a prefix that does not exist, you will not get any output.
    For example, this will not return anything since there is no project ID `SCPCP000099`:

    ```bash
    # List all files in a directory that does not exist
    aws s3 ls s3://openscpca-data-release/2024-05-01/SCPCP000099/
    ```

### How can I use results from existing modules in my analysis module?

Once modules reach maturity, the Data Lab adds them to a [Nextflow workflow](../ensuring-repro/openscpca-nf/index.md) to automatically generate the module's results from the latest data release.
If you want to use these module results (e.g., cell type annotations conducted in a different module) as input for your analysis, you can obtain them with [the `download-results.py` script](../getting-started/accessing-resources/getting-access-to-data.md#accessing-scpca-module-results).
Results downloaded with this script [will be stored in `data/current/results/{module name}`](../getting-started/accessing-resources/getting-access-to-data.md#downloaded-results-file-structure).

However, there may be circumstances when you want to use results from a module which has not yet been added to the Nextflow workflow.
In such cases, you will need to run the module yourself to generate the results.
Instructions for [running the module](../contributing-to-analyses/analysis-modules/running-a-module.md), including its software and compute requirements, should be available in the module's main `README.md` file.
After running the module, results will generally be stored in `analysis/{module name}/results`, and the module's documentation should describe the contents of result files.


### What if I want to use Seurat?

While [data downloads](../getting-started/accessing-resources/getting-access-to-data.md) are only available in `SingleCellExperiment` and `AnnData` format, `Seurat` versions of all objects (in [v5 assay format](https://satijalab.org/seurat/articles/seurat5_essential_commands)) are also available for use.

These files are part of the OpenScPCA results, associated with the module [`seurat-conversion`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/seurat-conversion) which we wrote to convert the processed `SingleCellExperiment` objects to `Seurat` format.
For more information on obtaining result files, please refer to the documentation for [the `download-results.py` script](../getting-started/accessing-resources/getting-access-to-data.md#accessing-scpca-module-results).

When working with these `Seurat` objects, please bear in mind the following:

* They were _not_ processed with a `Seurat` pipeline.
They were processed using the same pipeline as all OpenScPCA objects (e.g., with `Bioconductor`), and then converted to a `Seurat` format
  * Notably, they do contain the raw data counts, allowing you to perform normalization, dimension reduction, etc. with `Seurat` directly if you so choose
* To be more consistent with `Seurat` analysis pipelines, gene names in these objects use gene symbols rather than Ensembl ids


### The ScPCA data objects contain Ensembl ids, but I need gene symbols for my analysis. How should I perform this conversion?

In an effort to keep this consistent across the OpenScPCA project, we provide functions to convert Ensembl ids to gene symbols in an R package we maintain called [`rOpenScPCA`](https://github.com/AlexsLemonade/rOpenScPCA/).

This package has two particular functions to support this task:

* `rOpenScPCA::sce_to_symbols()`
  * This function converts row names in a `SingleCellExperiment` object from Ensembl ids to gene symbols
* `rOpenScPCA::ensembl_to_symbol()`
  * This function converts a vector of Ensembl ids to a vector of gene symbols

Please refer to these functions' help menus (e.g., `?rOpenScPCA::sce_to_symbols`) for additional information on how to use them.

### I noticed there are cluster assignments in the processed data files. Should I use those or re-cluster the data?

All ScPCA data objects contain cluster assignments which were [calculated using an automated pipeline](https://scpca.readthedocs.io/en/stable/processing_information.html#processed-gene-expression-data).
Because the clustering parameters used in this automated pipeline were not tailored to any given dataset, we do not recommend relying on these clusters for downstream analysis.
Instead, we strongly recommend re-clustering the data _and_ evaluating your cluster assignments before using them.

To support clustering analysis and evaluation, we provide several functions in an R package we maintain called [`rOpenScPCA`](https://github.com/AlexsLemonade/rOpenScPCA/) to accomplish the following tasks:

* Perform graph-based clustering
* Evaluate clustering results with quality control metrics
* Calculate several sets of clustering results across parameter space

We also provide an OpenScPCA analysis module [`hello-clusters`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-clusters) with example notebooks demonstrating how to use clustering functionality in `rOpenScPCA`.
