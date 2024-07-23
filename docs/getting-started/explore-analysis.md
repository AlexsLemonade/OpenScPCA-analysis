# Exploring existing analyses

Before contributing to OpenScPCA, you can explore and run any existing [analysis modules](../contributing-to-analyses/analysis-modules/index.md) on your local machine.

You can browse existing analysis modules in the repository in the [`analyses` folder](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses).
This folder also contains a pair of example analyses, one for [performing analysis in R](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-R) and one for [performing analysis in Python](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-python).

To run the example analysis modules or any other existing modules, you will need to:

- [Set up your local environment](#set-up-your-local-environment)
- [Follow instructions to run existing modules](#follow-instructions-to-run-existing-modules)
- [Get access to data for running modules](#get-access-to-data-for-running-modules)

## Set up your local environment

Before you can run any analysis modules locally, you need to set up your environment.

Please [refer to these linked instructions](../technical-setup/index.md) for the specific setup steps you need to take.

## Follow instructions to run existing modules

Each existing module has its own folder in the repository's [`analyses` folder](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses).
See our documentation on [analysis modules](../contributing-to-analyses/analysis-modules/index.md) to learn more about how analyses are structured.

All instructions and files ([except data](#get-access-to-data-for-running-modules)) needed to run a given module are provided in the module's folder.

- To [run existing modules](../contributing-to-analyses/analysis-modules/running-a-module.md), please refer to the `README.md` found in the root directory of that module.
- The `README.md` will contain a description of the module, instructions to run the module, the list of input and output files, and any additional software or computing resources you may need.


## Get access to data for running modules

For modules that use data obtained from the ScPCA Portal as input, you can directly download data from the [ScPCA Portal](https://scpca.alexslemonade.org/).

Additionally, you can use the [`simulate-sce`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/simulate-sce) analysis module to create simulated data:

- This module simulates data based on the properties of real ScPCA data.
- The simulated data contains permuted labels and simulated counts data for a small number of cells, making it ideal for testing analyses modules.
- The module requires either `SingleCellExperiment` or `AnnData` objects from the ScPCA Portal as input.

In some cases, an analysis module's input data will be the results/output from a different module.
These results are stored in an [S3 bucket on Amazon Web Services](../aws/index.md#s3-data-and-results-storage-with-aws), and are available to contributors.
Read more about [getting access to data](accessing-resources/getting-access-to-data.md).
