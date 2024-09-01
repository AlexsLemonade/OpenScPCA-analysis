This directory contains GitHub Actions (GHA) yaml files as well as helper files they use.

## Module testing

* Each workflow `run_{module name}.yml` is used to run the given module in CI using the test data.
* `run_all-modules.yml` runs all modules in CI by calling each `run_{module name}.yml` workflow.

## Building and pushing Docker images

* Each workflow `docker_{module name}.yml` is used to build the given module's Docker image.
* `docker_all-modules.yml` builds and pushes all module Docker images by calling each `docker_{module name}.yml` workflow.
* The `build-push-docker-module.yml` workflow contains actions used by `docker_all-modules.yml` and each `docker_{module name}.yml`.
It is not meant to be run as a standalone action.

## Creating releases

* `file_periodic_release_issues.yml` files issues to prepare for a periodic (dated) release.
* `create_periodic_release.yml` creates a periodic (dated) release.

## Repository maintenance

* `spellcheck.yml` performs spell checking of markdown and R Markdown documents.
* `code-styling.yml` performs code styling using `styler` and `ruff`.
  * This workflow also consumes `modules_code-styling.txt` to determine which modules to style.
  When modules reach the development stage where they should be automatically styled, they should be added to this text file.
