# Module workflows

To maintain module functionality over time, we use [GitHub Actions](https://docs.github.com/en/actions) (GHAs) to periodically run each module with the goal of testing that it runs to completion without errors.

For examples of existing analysis module GHAs, see the example Python and R module GHAs, [`run_hello-python.yml`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/.github/workflows/run_hello-python.yml) and [`run_hello-R.yml`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/.github/workflows/run_hello-R.yml), respectively.



All modules have a workflow which gets run upon changes to that module and periodically as cron. To begin it will use renv/conda to manage the environment but as the analysis matures, we'd like for it to use the docker file (note the timing of when to file such changes, link to docker section)

## Writing a module GHA

The Data Lab will generally step in to maintain and write these files, but you are also welcome to do so if you have the expertise!

When you [create a new module](../../contributing-to-analyses/analysis-modules/creating-a-module.md), an GHA workflow file is created in the file `.github/workflows/run_{module-name}.yml`.
Depending on the flags
On creation, this GHA does not yet run on any automatic triggers.

To prepare a module GHA for use, there are two



To make GHAs run efficiently, we strongly encourage and generally require, when appliable, that GHAs run the module code using the [simulated test data](../../getting-started/accessing-resources/getting-access-to-data.md#accessing-test-data).

This means that it's important to write your module code with sufficient flexibility to allow for test data to be used.
For example, when reading input data into scripts, you should use input arguments to specify input data paths, which will allow for test data input paths to be used in workflow modules.

