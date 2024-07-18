# Project tools

The OpenScPCA project supports UNIX-based platforms and is largely structured around the use of Git (and [GitHub](https://github.com)).
Therefore, as a contributor, there are several tools you will need to be familiar with, including:

- Git
    - [This documentation section](../../contributing-to-analyses/working-with-git/index.md) explains how to use Git as a contributor, specifically using the [GitKraken GUI](../../technical-setup/install-a-git-client.md#why-use-gitkraken)
- The [terminal](./using-the-terminal.md), also known as the command line
- Writing text in [markdown](./writing-in-markdown.md)
- Software management platforms, including [conda](../../ensuring-repro/managing-software/using-conda.md) and, if you plan to develop R-based modules, [`renv`](../../ensuring-repro/managing-software/using-renv.md)
    - We use conda and `renv` to manage installation of module-specific software dependencies and packages
- We additionally use [Docker](../../ensuring-repro/docker/index.md) and a [separate Nextflow workflow](../../ensuring-repro/openscpca-nf/index.md) to ensure reproducibility, but OpenScPCA contributors are not required to use these platforms unless they choose

!!! tip

    Are there other project tools you'd like to see documentation for?

    Let us know by [filing an issue to request additional documentation](https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/new?assignees=&labels=docs-request&projects=&template=04-docs-request.yml&title=Docs+request%3A)!
