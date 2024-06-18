# Automated module testing

To maintain module functionality over time, we use [GitHub Actions](https://docs.github.com/en/actions) (GHAs) to periodically run each module ("module testing GHA") with the goal of testing that the module code runs to completion without errors.

!!! info
    For more information about how we run modules to generate official results and OpenScPCA releases, please see our documentation on the `OpenScPCA-nf` workflow. <!-- openscpca-nf STUB_LINK -->

Module testing GHAs are automatically run in two circumstances:

- When a pull request is filed with changes to any module files
    - This GHA will need to pass without errors for [pull requests](../../contributing-to-analyses/pr-review-and-merge/index.md) to be approved
- On a periodic schedule
  - This ensures that changes in data or other code do not break tests within each module

For examples of existing analysis module GHAs, see the example Python and R module GHAs, [`run_hello-python.yml`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/.github/workflows/run_hello-python.yml) and [`run_hello-R.yml`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/.github/workflows/run_hello-R.yml), respectively.

To make GHAs run efficiently, the tests should run the module code with the [simulated test data](../../getting-started/accessing-resources/getting-access-to-data.md#accessing-test-data).
This means that it's important to write your module code with sufficient flexibility to allow for test data to be used.
You should read in files from the `data/current` directory, which will be automatically directed to test data during module testing GHA runs.

In addition, it's also helpful for your module to have a single entry point for running all module scripts and/or notebooks in their intended order, e.g. a [shell script](../../contributing-to-analyses/analysis-modules/running-a-module.md).
This way, the module testing GHA can directly call this script to execute the entire module.


## Writing a module testing GHA

!!! tip
    The Data Lab will generally maintain and write module testing GHAs, but you are welcome to do so as well if you are interested!
    See this GitHub documentation to learn about [workflow syntax for GHAs](https://docs.github.com/en/actions/using-workflows/workflow-syntax-for-github-actions).

When you [create a new module](../../contributing-to-analyses/analysis-modules/creating-a-module.md), a module testing GHA workflow file is created in the file `.github/workflows/run_{module-name}.yml`.
This initial file is inactive, meaning it will not run automatically run on the two aforementioned triggers.
As analysis module matures over time, the Data Lab staff will activate this GHA file so the module can be regularly tested.

### GHA steps

Each module testing GHA is initially created with these steps, which should be updated to reflect the given module's needs:

- Checkout the repository
- Download test data
    - Use the [`download-data.py`](../../getting-started/accessing-resources/getting-access-to-data.md#using-the-download-data-script) and/or [`download-results.py`](../../getting-started/accessing-resources/getting-access-to-data.md#accessing-scpca-module-results) scripts to specify the set of input files you need, with the `--test-data` flag to specify downloading the test data.
    - After this step, the `data/current` directory will point to the test data, ensuring the module GHA runs using the test data.
- Set up the module environment
    - Depending on [the flags used when creating your module](../../contributing-to-analyses/analysis-modules/creating-a-module.md#module-creation-script-flags), these steps will install the [`renv` and/or conda environment](../managing-software/index.md) from existing environment files (`renv.lock` and/or `conda-lock.yml`, respectively).
- Run the analysis module
    - Generally, this will involve calling the [module's run script](../../contributing-to-analyses/analysis-modules/running-a-module.md).

As an analysis module matures, the GHA will be updated to run the analysis in the [module's Docker image](../docker/docker-images.md), rather than using the `renv` and/or conda environment files.
Module testing GHAs can use their module's Docker images once the image has been built and pushed to the registry. <!-- STUB LINK building/updating docker images -->
