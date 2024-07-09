# OpenScPCA GitHub Action workflows

OpenScPCA uses [GitHub Actions](https://docs.github.com/en/actions) to automatically run several workflows that support overall reproducibility.

This section highlights two of these workflows:

1. [A module-testing workflow](./run-module-gha.md) which ensures that analysis modules can run without errors
2. [A Docker-building workflow](./build-docker-gha.md) which builds and pushes module-specific Docker images to a public registry

Generally, the Data Lab maintains and writes these GHA workflow files, but you are welcome to contribute as well if you are interested.

## Workflows in your fork

As described in each link above, workflow files are created when the [`create-analysis-module.py` script](../../contributing-to-analyses/analysis-modules/creating-a-module.md#module-workflows) is run.
The associated workflows themselves are inactive when they are first created.

Over the course of module development, these workflows will be activated so that they automatically run in several circumstances, including when pull requests are merged into the upstream `AlexsLemonade/OpenScPCA-analysis` repository's `main` branch.

GitHub will also run these activated workflows in your fork each time you [sync your fork's `main` branch](../../contributing-to-analyses/working-with-git/staying-in-sync-with-upstream.md) to the upstream `main` branch.
**When these workflows run in your fork, they will fail - this is expected and not a concern!**
The reason they fail is because they use credentials that are, by design, only available to the `AlexsLemonade/OpenScPCA-analysis` repository.

If you prefer to turn off GitHub Actions in your fork to avoid these failures, you can update your fork's settings:

- In the repository Settings menu, select `Actions -> General` from the left-hand menu
- Turn on the setting "Disable actions"
