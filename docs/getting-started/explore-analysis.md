# Explore existing analysis

Before contributing to OpenScPCA, you can explore and run any existing [analysis modules](../contributing-to-analyses/analysis-modules/index.md) on your local machine.

You can browse existing analysis modules in the repository in the [`analyses` folder](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses).
This folder also contains a pair of example analyses, one for [performing analysis in R](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-R) and one for [performing analysis in Python](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/hello-python).

To run the example analysis or any other existing analyses, you will need to:

- [Set up local environment](#set-up-local-environment)
- [Run existing modules](#run-existing-modules)
- [Accessing data for modules](#accessing-data-for-modules)

## Set up local environment

Before you can run any analysis modules locally, you need to [set up your environment](../technical-setup/index.md).
This includes:

- Downloading and setting up a [Git client](../technical-setup/install-a-git-client.md)
- [Forking the `AlexsLemonade/OpenScPCA-analysis` repository](../technical-setup/fork-the-repo.md)-
- [Cloning your fork](../technical-setup/clone-the-repo.md) to your computer
- [Installing conda and necessary dependencies](../technical-setup/environment-setup/index.md) needed to run analysis

## Run existing modules

Each existing module will have it's own folder in `analyses`.
All instructions and files needed to run a given module will be contained in the module's folder.

- To run existing modules, please reference the `README.md` found in the root directory for that module.
- The `README.md` will contain a description of the module, instructions to run the module, the list of input and output files, and any additional software or computing resources you may need.

See [analysis module](../contributing-to-analyses/analysis-modules/index.md) to learn more about how analyses are structured.

## Accessing data for modules

The specific requirements for input data will be listed in the `README.md` for each module.
For modules that use data obtained from the ScPCA Portal as input, data can be directly downloaded from the [ScPCA Portal](https://scpca.alexslemonade.org/).

Additionally, you can use the [`simulate-sce`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/simulate-sce) analysis module to create simulated data.

- This module simulates data based on the properties of real ScPCA data.
- The simulated data contains permuted labels and simulated counts data for a small number of cells, making it ideal for testing analyses modules.
- The module requires either `SingleCellExperiment` or `AnnData` objects from the ScPCA Portal as input.

Any modules that use input data obtained from another module are stored in an S3 bucket on Amazon Web Services and are available to contributors.
<!--TODO: Fill in with link to getting access to data-->
More information on getting access to data coming soon!
