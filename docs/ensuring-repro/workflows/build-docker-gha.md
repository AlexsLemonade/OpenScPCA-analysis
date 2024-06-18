# Automated Dockerfile building

We use [GitHub Actions](https://docs.github.com/en/actions) (GHAs) to build Dockerfiles and push images to [Amazon ECR](https://aws.amazon.com/ecr/), a registry of containers such as Docker images.

Docker building GHAs are automatically run in two circumstances:

- When a pull request is filed with changes to any module's environment files, such as `renv.lock`, `conda.lock`, or the `Dockerfile` itself
- When new `OpenScPCA-analysis` releases are made


These pushed images are used in a variety of ways:

- The `OpenScPCA-nf` workflow pulls module-specific Docker images to reproducibly run modules and generate results <!-- STUB_LINK openscpca-nf -->
- [Module testing GHAs](./run-module-gha.md) will, as the given analysis module matures, pull the Docker-specific Docker image to create the environment used in the workflow
- OpenScPCA contributors, as well as the wider research community, can freely pull module-specific images to reproducibly run OpenScPCA analysis modules, for example to locally run analysis modules or to develop within a Docker container

For examples of existing Docker building GHAs, see the example `simulate-sce` GHA [`docker_simulate-sce.yml`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/.github/workflows/docker_simulate-sce.yml).


## Writing a Docker building GHA

!!! tip
    The Data Lab will generally maintain and write Docker building GHAs, but you are welcome to do so as well if you are interested!
    See this GitHub documentation to learn about [workflow syntax for GHAs](https://docs.github.com/en/actions/using-workflows/workflow-syntax-for-github-actions).

When you [create a new module](../../contributing-to-analyses/analysis-modules/creating-a-module.md), a Docker building GHA workflow file is created in the file `.github/workflows/docker_{module-name}.yml`.
This initial file is inactive, meaning it will not run automatically run on the two aforementioned triggers.
As analysis module matures over time, the Data Lab staff will activate this GHA file so the Docker image can be built and pushed to the Amazon ECR registry for general use.
