# OpenScPCA-analysis

OpenScPCA is an open, collaborative project to analyze data from the [Single-cell Pediatric Cancer Atlas (ScPCA) Portal](https://scpca.alexslemonade.org/), which currently holds 500 samples from over 50 pediatric cancer types.

This project aims to:

- Characterize the ScPCA data with analyses such as labeling cell types or identifying recurrent cell states in multiple tumor types

- Work on open and collaborative analyses

- Build consensus around usage, strengths, and pitfalls of methods and their application to pediatric cancer data

- Improve the utility of the ScPCA data for the research community by generating new data representations (e.g., integrated objects)

### Platform and Language

The primary platform for the project is Linux, and we're best able to support analyses written in R and Python.
The documentation and analyses modules cater to Linux and Mac platforms, and the above languages.

If desired, please [file an issue](https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/new?assignees=&labels=docs-request&projects=&template=04-docs-request.yml&title=Docs+request%3A) to request Windows support and associated docs.

## Contributing

To start contributing:

1. Please review our [Policies](https://openscpca.readthedocs.io/en/latest/policies/).

2. Fill out the [interest form](https://share.hsforms.com/1MlLtkGYSQa6j23HY_0fKaw336z0).

3. Visit [Getting Started](https://openscpca.readthedocs.io/en/latest/getting-started/making-your-first-analysis-contribution/) for first steps.

### Communicating

You can ask questions, propose analyses, get help in [GitHub Discussions](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions).

You can also join the OpenScPCA community in the `#open-scpca` channel on our [Cancer Data Science Slack](https://ccdatalab.org/slack).

GitHub Issues are reserved for planned and defined tasks or bug reports.
If this is your first time interacting with the project, please post in Discussions.

If you wish to report a security vulnerability, please [email](mailto:report@ccdatalab.org) us.
Do NOT report it in a public forum.
See our [security policy](./SECURITY.md) for more information.

Please see our documentation on [Tools for communication](https://openscpca.readthedocs.io/en/latest/communications-tools/) for more information.

## Documentation

We have comprehensive documentation to help you with various aspects of participating in OpenScPCA - from proposing an analysis, acquiring data, troubleshooting errors to making your first code contribution.

Please refer to it as you work on the project.

|[OpenScPCA Documentation](https://openscpca.readthedocs.io/en/latest/)|
|---|

## Setting up and running analyses

### Local Setup

To set up, you will need a Git client, Miniconda, R or Python, and an editor of your choice. Please see our documentation on [Technical Setup](https://openscpca.readthedocs.io/en/latest/technical-setup) for detailed instructions to create a local setup.

### Accessing Data

Data is publicly available from the [ScPCA Portal](https://scpca.alexslemonade.org/) and from an AWS S3 bucket for OpenScPCA project contributors.
We also provide a way to download smaller, simulated data files for you to play with.

Please refer to [Getting Access to Data](https://openscpca.readthedocs.io/en/latest/getting-started/accessing-resources/getting-access-to-data/) for more details.

### Running an analysis module

Each analysis module has a `README.md` file which contains instructions to run that specific module.
Please see the relevant analysis module's `README.md` for instructions.

Please see our documentation on [running an analysis module](https://openscpca.readthedocs.io/en/latest/contributing-to-analyses/analysis-modules/running-a-module/) for more information.
