# Using Docker images

Once a module's Docker image has been [pushed to the OpenScPCA Docker registry in the Amazon ECR Public Gallery](#STUB_LINK/workflows/build-docker-gha.md), you can obtain the image to locally run modules in a reproducible environment.
While modules can also be run within specified [conda and/or `renv` environments](../managing-software/index.md), using Docker images provides a more streamlined (but advanced) approach to ensuring a fully reproducible software environment.

The OpenScPCA Docker register is publicly available at <https://gallery.ecr.aws/openscpca/> and contains Docker images with tags of the form `public.ecr.aws/openscpca/{module-name}:latest`.

!!! note
      This page presents the main steps for obtaining and running a module in a Docker image.
      These instructions assume that you have already [downloaded Docker Desktop](./index.md#how-to-install-docker) and are familiar with using Docker images and launching containers.

      To learn more about using Docker, see our recommended resources [on this page](./index.md).

      If you would like to see additional dociumentation about working with Docker images, please [let us know by filing an issue](https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/new?assignees=&labels=docs-request&projects=&template=04-docs-request.yml&title=Docs+request%3A).

## Obtaining Docker images

You can see all available OpenScPCA Docker images at this URL: <https://gallery.ecr.aws/openscpca/>.
Images are named `openscpca/{module-name}`.

To obtain a local copy of a Docker image, follow these steps:

1. Open the Docker Desktop application and ensure it is running.
2. From a [terminal](../../getting-started/project-tools/using-the-terminal.md), pull the image using this command, where `{module-name}` is replaced with the given Docker image you wish to pull:

      ```sh
      docker pull public.ecr.aws/openscpca/{module-name}:latest
      ```

    !!! tip "Macs with Apple silicon"
          If you are using an Apple silicon (M1-3) Mac, you will need to use the additional flag `--platform linux/x86_64` when pulling an image:

          ```sh
          docker pull --platform linux/x86_64 public.ecr.aws/openscpca/{module-name}:latest
          ```

3. We recommend updating your Docker Desktop resource settings to ensure Docker has access to sufficient compute resources to run the module. For more information on module compute requirements, please refer to the [module's `README.md` file](../../contributing-to-analyses/analysis-modules/compute-requirements.md#readme-files).
      - Follow these instructions to update your Docker Desktop settings for [Mac](https://docs.docker.com/desktop/settings/mac#resources) and [Windows](https://docs.docker.com/desktop/settings/windows/#resources) computers. You will need to click "Apply and Restart" to ensure these changes go through.

4. You can now use the image to run code from the analysis module using [the `docker run` command](https://docs.docker.com/reference/cli/docker/container/run/) from a terminal.
      - The analysis module's `README.md` file should contain documentation for what command(s) you need to issue to [run the module](../../contributing-to-analyses/analysis-modules/running-a-module.md).
      - You will also need to [mount a volume](https://docs.docker.com/storage/volumes/) as part of your command to ensure that the Docker container has access to all files needed, including any data files.
      - Note that containers for [images based on `bioconductor/r-ver`](./docker-images.md#r-based-images) can also be run in the browser as interactive RStudio Server sessions, [as described in the Bioconductor documentation](https://www.bioconductor.org/help/docker/).


