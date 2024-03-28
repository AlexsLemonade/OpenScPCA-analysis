# OpenScPCA-analysis

OpenScPCA is an open, collaborative project to analyze data from the [Single-cell Pediatric Cancer Atlas (ScPCA) Portal](https://scpca.alexslemonade.org/), which currently holds 500 samples from over 50 pediatric cancer types. 

This project aims to: 

- Characterize the ScPCA data with analyses such as labeling cell types or identifying recurrent cell states in multiple tumor types

- Work on open and collaborative analyses.

- Build consensus around usage, strengths, and pitfalls of methods and their application to pediatric cancer data.

- Improve the utility of the ScPCA data for the research community by exploring the creation of new data assets (e.g., integrated objects).

### Platform and Language
The primary platform for the project is Linux and we're best able to support analyses written in R and Python.
The documentation and analyses modules cater mainly to Linux and Mac platforms, and the above languages.

If you would like support for Windows, please [file an issue](https://github.com/AlexsLemonade/OpenScPCA-analysis/issues/new?assignees=&labels=docs-request&projects=&template=04-docs-request.yml&title=Docs+request%3A) to request Windows support and associated docs.

## Contributing

To start contributing:

1. Please review our [policies](#STUB_LINK)

2. Fill out the interest form (TODO: make case for interest form)

3. Visit [Getting Started](#STUB_LINK) for first steps. 

### Communicating

You can ask questions, propose analyses, get help in [Github Discussions](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions).

You can also join the OpenScPCA community in the #open_scpca channel on our [Cancer Data Science](slack).

Please see [Tool for communication](#STUB_LINK) for more information.

## Setting up and running analyses

### Local Setup

To setup, you will need a Git client, mini conda, R or Python and a editor of your choice. Please see [Technical Setup] for detailed instructions to create a local setup.

### Accessing Data

Data is available via the [ScPCA Portal](https://scpca.alexslemonade.org/) and also via an AWS S3 bucket.
We have also provided a way to generate simulated data for you to play with. 

Please [Getting access to Data](#STUB_LINK) for more details.

### Running a module

Each module has a README which contains instructions to run that specific module. Please see the relevant analysis module's README for instructions.