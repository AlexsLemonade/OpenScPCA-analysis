# Automated Dockerfile building and pushing

We use [GitHub Actions](https://docs.github.com/en/actions) (GHAs) to build Dockerfiles and push images to [Amazon ECR](https://aws.amazon.com/ecr/), a registry of pre-built Docker images.

We generally build and/or push Docker images to ECR under two circumstances:

1. When a pull request is filed with changes to the given module's environment files, such as `renv.lock`, `conda.lock`, or the `Dockerfile` itself
    - The pull request itself triggers a GHA to test that the Dockerfile builds
    - When the pull request is merged, the Docker image is again built and pushed to ECR
2. When a new `OpenScPCA-analysis` release is made
    - `OpenScPCA-analysis` releases automatically trigger building and pushing all module-specific Docker images to ECR

After they have been pushed, the Docker images are used in a variety of ways:

- [Module testing GHAs](./run-module-gha.md) may pull the module-specific Docker image to run the module code in the pre-built image environment
- The `OpenScPCA-nf` workflow pulls module-specific Docker images to reproducibly run modules and generate results <!-- STUB_LINK openscpca-nf -->
- OpenScPCA contributors, as well as the wider research community, can freely pull module-specific images to reproducibly run OpenScPCA analysis modules, for example to locally run analysis modules

All images pushed to ECR will be available from the [OpenScPCA Docker registry in the Amazon ECR Public Gallery](https://gallery.ecr.aws/openscpca) with image tags of the form `public.ecr.aws/openscpca/{module-name}:latest`.
For more information about pulling Docker images from ECR and using them locally, please see our documentation on [using Docker images](../docker/using-images.md).

## Writing a Docker building GHA

!!! tip
    The Data Lab will generally maintain and write Docker building GHAs, but you are welcome to do so as well if you are interested!
    See this GitHub documentation to learn about [workflow syntax for GHAs](https://docs.github.com/en/actions/using-workflows/workflow-syntax-for-github-actions).

When you [create a new module](../../contributing-to-analyses/analysis-modules/creating-a-module.md), a Docker building GHA workflow file is created in the file `.github/workflows/docker_{module-name}.yml`.
This initial file is inactive, meaning it will not run automatically run on the two aforementioned triggers.
Once a [Dockerfile with module dependencies has been added to a given analysis module](../docker/docker-images.md#analysis-module-dockerfiles), the Data Lab staff will [activate this GHA file](../docker/docker-images.md#automated-testing-and-deployment) so the Docker image can be built and pushed to the Amazon ECR registry for general use.
