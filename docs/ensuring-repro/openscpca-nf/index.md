# The OpenScPCA Nextflow workflow

In addition to the `OpenScPCA-analysis` repository, the Data Lab maintains a separate [`OpenScPCA-nf` repository](https://github.com/AlexsLemonade/OpenScPCA-nf) that holds a [Nextflow workflow](https://www.nextflow.io/) which we use to generate results from mature `OpenScPCA-analysis` modules.

As analysis modules are completed, the Data Lab will port the analysis code into the `OpenScPCA-nf` repository and add it to the Nextflow workflow.
To ensure reproducibility, the workflow will run each module [in its Docker image](../docker/docker-images.md).

## `OpenScPCA-nf` workflow results

The Data Lab will run the `OpenScPCA-nf` workflow on the same schedule as we make `OpenScPCA` data releases, and/or if there have been sufficient changes to the `OpenScPCA-nf` workflow to merit re-generating results.
We will use the workflow to generate two sets of results for each completed module:

- Results generated from the real `OpenScPCA` data
- Results generated from the [simulated test `OpenScPCA` data](../../getting-started/accessing-resources/getting-access-to-data.md#accessing-test-data)

In some cases, you may wish to use results from a completed module in your analysis.
For example, you may wish to use cell type labels from a module that performed cell type annotation on a set of ScPCA samples.

You can obtain either of these `OpenScPCA` results using the provided [`download-results.py` script](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/download-results.py).
Please refer to our [documentation on the `download-results.py` script](../../getting-started/accessing-resources/getting-access-to-data.md#accessing-scpca-module-results) for more information on its usage.

Note that if you need to use results from a module has not yet been added to the `OpenScPCA-nf` workflow, you will need to [run the module yourself](../../contributing-to-analyses/analysis-modules/running-a-module.md) to generate results.
