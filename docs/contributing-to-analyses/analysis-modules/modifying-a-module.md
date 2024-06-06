# Modifying an existing module


If you are modifying a module that already exists, then you will first need to file an issue using the [`Update an analysis` issue template](../../communications-tools/github-issues/issue-templates.md#update-an-analysis).
Then you can create a new branch and commit any code changes before filing a PR.

Be sure that you add documentation for any code changes that you make, as well as any new software dependencies your changes require.
You will also need to [upload the updated module results to S3](../../software-platforms/aws/working-with-s3-buckets.md#syncing-your-results-to-s3).

You should include these documentation changes in any PRs that modify existing modules:

- Fully document your changes [in the module's `README.md` file](./documenting-analysis.md).
This includes all of:
    - Changes to the behavior of any existing script or notebooks you modify
    - New scripts or notebooks you add to the module
    - Changes to which input data the module consumes
    - Changes to any output files (e.g., results or plots) the module generates
    - Changes to the module's [compute requirements](./compute-requirements.md)
    - Changes to instructions for setting up the module's run environment and/or running the module
- Make sure any new dependencies have been added to the module's [software environment](../../ensuring-repro/managing-software/index.md).
- Either update or add `README.md` files to any sub-folders you modify or add to the module.
    - For example, if you added a new script to the `scripts` directory, be sure to add a `scripts/README.md` that describes the usage of the new script.
