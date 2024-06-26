# Using Docker images

Each OpenScPCA analysis module has [its own Docker image](./docker-images.md) to create a fully reproducible environment for running the module.
Once a module's Docker image has been [pushed to the OpenScPCA Docker registry in the Amazon ECR Public Gallery](../workflows/build-docker-gha.md), you can obtain the image to locally run modules in a reproducible environment.
While modules can also be run within specified [conda and/or `renv` environments](../managing-software/index.md), using Docker images provides a more streamlined (but advanced) approach to ensuring a fully reproducible software environment.

!!! note
      This page presents the main steps for obtaining and running a module in a Docker image.
      These instructions assume that you have already [downloaded Docker Desktop](./index.md#how-to-install-docker) and are familiar with using Docker images and launching containers.

      To learn more about using Docker, see our recommended resources [on this page](./index.md).

      If you would like to see additional dociumentation about working with Docker images, please [let us know by filing an issue](https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/new?assignees=&labels=docs-request&projects=&template=04-docs-request.yml&title=Docs+request%3A).

## Obtaining Docker images

The OpenScPCA Docker register is publicly available at <https://gallery.ecr.aws/openscpca/> and lists all available OpenScPCA Docker images, which are named with tags of the form `public.ecr.aws/openscpca/{module-name}:latest`.

!!! tip "All Docker images have lowercase names"
      OpenScPCA Docker images are named `openscpca/{module-name}`, where `{module-name}` is a lowercase version of the module name.
      For example, the Docker image for a module called `hello-R` would be named `openscpca/hello-r`.

To obtain a local copy of a Docker image, follow these steps:

1. Open the Docker Desktop application and ensure it is running.
2. From a [terminal](../../getting-started/project-tools/using-the-terminal.md), pull the image using this command, where `{module-name}` is replaced with the given Docker image you wish to pull:

      ```sh
      docker pull public.ecr.aws/openscpca/{module-name}:latest
      ```

    !!! tip "Macs with Apple silicon"
          If you are using an Apple silicon (M1-3) Mac, you will need to use the additional flag `--platform linux/x86_64` when pulling an image.
          Note that the flag `--platform linux/x86_64` _must_ be provided before the image name.

          ```sh
          docker pull --platform linux/x86_64 public.ecr.aws/openscpca/{module-name}:latest
          ```

## Running Docker images

### Update your resource settings

Before running any Docker image, we recommend updating your Docker Desktop resource settings to ensure Docker has access to sufficient compute resources to run the module.
For more information on module compute requirements, please refer to the [module's `README.md` file](../../contributing-to-analyses/analysis-modules/compute-requirements.md#readme-files).

Follow these instructions to update your Docker Desktop settings for [Mac](https://docs.docker.com/desktop/settings/mac#resources) and [Windows](https://docs.docker.com/desktop/settings/windows/#resources) computers.
You will need to click "Apply and Restart" to ensure these changes go through.

### Using `docker run`

You can now use the image to run the analysis module in a container using [the `docker run` command](https://docs.docker.com/reference/cli/docker/container/run/) from a terminal.
The analysis module's `README.md` file should contain documentation for what command(s) you need to issue to [run the module](../../contributing-to-analyses/analysis-modules/running-a-module.md).

When launching a container, you should mount the `OpenScPCA-analysis` repository directory as a volume in the Docker image, rather than just the analysis module directory itself.
This way, the Docker container will have access to the full project repository including the `data` directory where most module input data is stored.

For example, you might launch a Bash session in a container with the following line, which will mount your `OpenScPCA-analysis` directory into the container's `/home` directory.
Once launched, you can navigate to the analysis module of interest and run its code.

```sh
docker run \
  -it \
  --rm \
  --mount type=bind,source=/local/path/to/OpenScPCA-analysis,target=/home/OpenScPCA-analysis \
  public.ecr.aws/openscpca/{module-name}:latest \
  bash
```

!!! tip "Macs with Apple silicon"
      If you are using an Apple silicon (M1-3) Mac, you may see this warning when launching the container:

      ```{.console .no-copy}
      WARNING: The requested image's platform (linux/amd64) does not match the detected host platform (linux/arm64/v8) and no specific platform was requested
      ```

      You can safely ignore this warning, or you can silence it by providing the flag `--platform linux/x86_64` to your `docker run` command.


