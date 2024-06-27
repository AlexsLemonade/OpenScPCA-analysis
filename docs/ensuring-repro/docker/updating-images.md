# Updating Docker images

Over the course of module development, new dependencies and packages may get added or removed, thereby changing the [module's software environment](../managing-software/index.md).
For example, your environment will change if you add a new R or Python package to perform additional analyses.
Changing the software environment means that module-specific Docker image also needs to be updated to incorporate those changes.

!!! note
    These instructions assume that the Data Lab has already activated both of the [GitHub Action (GHA) workflows](../../contributing-to-analyses/analysis-modules/creating-a-module.md#module-workflows) that were created when the module was established, including the [module-testing GHA](../workflows/run-module-gha.md) and the [Docker-building GHA](#STUB_LINK../workflows/docker-build-gha.md).

    In addition, the module-testing GHA should have already been updated to build its runtime environment from the module's Docker image, not from other environment files (e.g. `renv.lock` or `conda-lock.yml`.)

    If these conditions are not met, you do not need to take any specific steps to update Docker images when you update a module's software environment.


## Context (header help????)

Docker images in OpenScPCA are used in two primary circumstances:

1. When generating module results as part of the `OpenScPCA-nf` workflow <!-- STUB_LINK -->
2. When the [module-testing GitHub Action (GHA) workflow](../workflows/run-module-gha.md) is run on pull requests (PRs) that modify module code
    - This GHA pulls the module-specific Docker image from the [OpenScPCA Docker registry](https://gallery.ecr.aws/openscpca/) to create a reproducible runtime environment with all module dependencies

The module-testing GHA will fail if any module dependencies are missing from the Docker image it pulls down.
Therefore, for this GHA to pass, any environment changes must have already been incorporated into the module's Docker image and pushed to the Docker registry.

## Steps to update a Docker image

Whenever you add or remove a new package or dependency to your module's software environment, you will generally need to take a two-step process to add new code that uses those new dependencies:

1. First, file a PR that _only contains_ your updated environment files, e.g. `renv.lock` and/or `conda-lock.yml`.
_This PR should not any additional code changes that need this updated environment._
    - Once this PR is merged, the module's Docker image will be [rebuilt and pushed to the OpenScPCA Docker Registry](#STUB_LINK../workflows/docker-build-gha.md)
1. Then, you can file one or more PRs with code changes that use the update software environment.
    - By ensuring the updated Docker image was pushed to the registry _before_ filing this PR, the [module testing GHA](../workflows/run-module-gha.md) will be able to use the most up-to-date Docker image and avoid dependency errors
