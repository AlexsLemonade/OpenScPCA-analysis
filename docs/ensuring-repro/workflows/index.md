# Automated OpenScPCA workflows

OpenScPCA uses [GitHub Actions](https://docs.github.com/en/actions) to automatically run several workflows that support overall reproducibility.

This section highlights two of these workflows:

1. [A module-testing workflow](./run-module-gha.md) which ensures that analysis modules can run without errors
2. [A Docker-building workflow](./build-docker-gha.md) which builds and pushes module-specific Docker images to a public registry

Generally, the Data Lab maintains and writes these GHA workflow files, but you are welcome to contribute as well if you are interested.
