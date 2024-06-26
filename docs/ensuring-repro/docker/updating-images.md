# Updating Docker images

Over the course of module development, new dependencies may get added or removed, thereby changing the [module's software environment](../managing-software/index.md).
As modules progress to maturity, there are important considerations to bear in mind when updating the module's software environment.

In particular, over the course of module development, the Data Lab will activate two [GitHub Action (GHA) workflows](../workflows/index.md) to ensure module reproducibility:

- A GHA workflow to [build and push the module's Docker image](#STUB_LINK../workflows/docker-build-gha.md) to the OpenScPCA Docker Registry
    - This GHA automatically runs when pull requests (PRs) that modify a module's environment are merged into `main`
- A GHA workflow to [run the module using test data](../workflows/run-module-gha.md)
    - This GHA automatically runs on PRs that modify any module code
    - Once the module's Docker image has been pushed to the OpenScPCA Docker Registry, this GHA will pull that Docker image to build a runtime environment.
    Therefore, for this GHA to succeed, it must use the most up-to-date Docker image.

Once both of these workflows have been activated, you will need to file two PRs every time you write code that takes advantage of changes to the module software environment:

1. First, file a PR that contains _only_ changes to your module's environment files, e.g. `renv.lock` and/or `conda.lock`.
This PR should not include additional code changes that require the updated environment.
    - Once this PR is merged, the module's Docker image will be rebuilt and pushed to the OpenScPCA Docker Registry
1. Second, file a PR that includes the code changes which use the update software environment
    - By ensuring the updated Docker image is pushed to the registry _before_ filing this PR, the [module testing GHA](../workflows/run-module-gha.md) will be able to use most up-to-date Docker image and avoid dependency errors
