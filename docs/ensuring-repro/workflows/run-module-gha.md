# Automated module testing

To maintain module functionality over time, we use [GitHub Actions](https://docs.github.com/en/actions) (GHAs) to periodically run each module ("module testing GHA") with the goal of testing that it runs to completion without errors.

!!! info
    For more information about how we run modules to generate official results and OpenScPCA releases, please see our documentation on the `OpenScPCA-nf` workflow. <!-- openscpca-nf STUB_LINK -->

Module testing GHAs are automatically run under these two circumstances:

- When a pull request is filed with changes to any module files
    - This GHA will need to pass without errors for [pull requests](../../contributing-to-analyses/pr-review-and-merge/index.md) to be approved
- All modules are periodically run on a schedule to ensure they continue to pass tests as OpenScPCA grows

To make GHAs run efficiently, we strongly encourage and generally require, when applicable, that GHAs run the module code using the [simulated test data](../../getting-started/accessing-resources/getting-access-to-data.md#accessing-test-data).
This means that it's important to write your module code with sufficient flexibility to allow for test data to be used.
For example, when reading input data into scripts, you should use arguments for input data files, allowing for test data paths to be used in these module testing GHAs.

For examples of existing analysis module GHAs, see the example Python and R module GHAs, [`run_hello-python.yml`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/.github/workflows/run_hello-python.yml) and [`run_hello-R.yml`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/.github/workflows/run_hello-R.yml), respectively.

## Writing a module testing GHA

!!! tip
    The Data Lab will generally maintain and write module testing GHAs, but you are welcome to do so as well if you are interested!
    See this GitHub documentation to learn about [workflow syntax for GHAs](https://docs.github.com/en/actions/using-workflows/workflow-syntax-for-github-actions).

When you [create a new module](../../contributing-to-analyses/analysis-modules/creating-a-module.md), a GHA workflow file is created in the file `.github/workflows/run_{module-name}.yml`.
This initial file is inactive, meaning it will not run automatically run on the two aforementioned triggers.
As analysis module begins to mature over time, the Data Lab staff will activate this workflow file so the module can be regularly tested.

### GHA steps

Each module testing GHA is initially created with these steps, which should be updated to reflect the given module's needs:

- Checkout the repository
- Download test data
- Setup the module environment
    - Depending on [the flags used when creating your module](../../contributing-to-analyses/analysis-modules/creating-a-module.md#module-creation-script-flags), this will steps steps needed to install the [`renv` and/or conda environment](../managing-software/index.md) from existing environment files (`renv.lock` and/or `conda-lock.yml`, respectively)
- Run the analysis module
    - Generally, this will involve calling the [module's run script](../../contributing-to-analyses/analysis-modules/running-a-module.md)

As an analysis module matures, the GHA should instead be updated to run the analysis in the module's Docker image, rather than using the `renv` and/or conda environment files.
Module testing GHAs can use their module's Docker images once the image has been built and pushed to the registry. <!-- STUB LINK building/updating docker images -->
