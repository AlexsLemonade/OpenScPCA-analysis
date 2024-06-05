# Getting access to data

To participate in the OpenScPCA Project, you will need access to the Single-cell Pediatric Cancer Atlas (ScPCA) data.
Here, we describe several ways to access data.

Broadly speaking, there are three kinds of data you might wish to work with:

- Data from the ScPCA Portal
    - You can find out more about the contents of files and how they were processed from the ScPCA documentation: <https://scpca.readthedocs.io>
- Simulated data for testing and/or developing modules
- Results from completed OpenScPCA analysis modules

Below, we describe how you can obtain each of these types of data.


## Accessing ScPCA Data

### Accessing data on the ScPCA Portal

ScPCA data are readily available from the ScPCA Portal that the Data Lab maintains: <https://scpca.alexslemonade.org/>.

You can select the project(s) or sample(s) you are interested in analyzing and download them from the Portal.

We recommend creating a `portal_downloads` subdirectory in the local copy of the `data` directory, which can be accomplished by running the following command from the root of the repository on Linux or macOS:

```sh
mkdir -p data/portal_downloads
```

You can then develop your analysis using these paths temporarily.

Before filing a pull request, you should change your paths to reflect the directory established by the download script (`data/current`) described in the next section.

### Accessing data from S3

Because we expect that contributors may want to work with many samples or analyze data on remote systems, ScPCA data is also available to contributors on AWS S3 via a script for downloading data.

Before you can access data in this manner, the Data Lab team needs to create an AWS account for you; these data are not publicly accessible.
See our documentation [getting access to AWS](index.md#getting-access-to-aws) for more information.

!!! info Review additional restrictions
    Please review the `DATA_USAGE.md` file included in all downloads for any additional restrictions on data usage (see [Policies](../../policies/index.md)).

#### Using the download data script

!!! note
    These instructions assume you have taken the following steps:

    - [Cloned the repo](../../technical-setup/clone-the-repo.md)
    - [Set up conda](../../technical-setup/environment-setup/setup-conda.md) (note that conda is pre-installed, but not yet set up, on [Lightsail for Research](../../software-platforms/aws/starting-development-on-lsfr.md#create-and-activate-a-conda-environment) instances)
    - [Configured the AWS CLI and logged in](../../technical-setup/environment-setup/configure-aws-cli.md) OR are [using Lightsail for Research](../../software-platforms/aws/index.md#lsfr-virtual-computing-with-aws) (which doesn't require logging in via the AWS CLI)

The [`download-data.py` script](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/download-data.py) is designed to download files from whatever release you specify to a folder named for the date of that release and [symlink](https://en.wikipedia.org/wiki/Symbolic_link) it to `data/current`.

We briefly cover some of the most common use cases below, but we encourage you to review all the options available to you.

You can list all the options for the download data script by running the following from the root directory of the repository:

```sh
./download-data.py --help
```

!!! tip "No `--profile` necessary on Lightsail for Research"
    Omit `--profile openscpca` from the commands below if you're using Lightsail for Research.

Assuming you created a profile called `openscpca` when [configuring the AWS CLI](../../technical-setup/environment-setup/configure-aws-cli.md), you can run the download script with all default options to download all processed samples from the most recent release in `SingleCellExperiment` format:

```sh
./download-data.py --profile openscpca
```

If you prefer to work with `AnnData` files, you can use the `--format` option as follows:

```sh
./download-data.py --format AnnData --profile openscpca
```

To review what samples would be downloaded without performing the download yet, you can use the `--dryrun` option as follows:

```sh
./download-data.py --dryrun --profile openscpca
```

If you're only working with a subset of the data, you can use the `--projects` or `--samples` to download select samples (note: these options are mutually exclusive).

#### Download data structure

Downloads via the download script will generally have the following structure:

```sh
data
├── {Release}
│       ├── {Project ID}
│       │   └── {Sample ID}
│       │       └── {Library files}
│       ├── bulk_metadata.tsv (if applicable)
│       ├── bulk_quant.tsv (if applicable)
│       └── single_cell_metadata.tsv
└── current -> {Release}
```

See the ScPCA documentation for [more information about individual files](https://scpca.readthedocs.io/en/latest/sce_file_contents.html).

Data releases are dated using the following format: `YYYY-MM-DD`.
By default, the data download script will download the most recent release.

To list all available releases, you can use the following command from the root of the repository:

```sh
./download-data.py --list-releases
```

The release directory is symlinked to `data/current`.
This is generally the path you should use in your code.

## Accessing simulated test data

You can also use the [download data script](#using-the-download-data-script) to download the test data for the project by using the `--test-data` option.
The test data are generally smaller, simulated or permuted data files with the same structure as the real project data (learn more at the [`simulate-sce`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/simulate-sce) module).

You do not need an AWS account set up to download the test data.

To download test data with all other options set at default, run the following from the root of the repository:

```sh
./download-data.py --test-data
```

This will download test data and direct `data/current` symlink to the test data directory.

To switch back to using ScPCA data, rerun the script with the `--release current` option.
If you already have the most recent data, this will not repeat downloading the data you already have.


## Accessing ScPCA module results

Once an analysis modules matures, the Data Lab ports the analysis to a separate repository called [`OpenScPCA-nf`](https://github.com/AlexsLemonade/OpenScPCA-nf).
This repository holds a [Nextflow workflow](https://www.nextflow.io/) to reproducibily run and generate results from completed analysis modules.

If you wish to use output from a completed analysis module in your module, you will need to obtain the results generated by this workflow.

We have provided a [`download-results.py` script](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/download-results.py) that you can use to download all results from a module that has been run in `OpenScPCA-nf`.
Results are available at specific tagged releases. <!-- STUB_LINK to openscpca-nf docs -->

We briefly cover some of the most common use cases below, but we encourage you to review all the options available to you.

You can list all the options for `download-results.py` by running the following from the root directory of the repository:

```sh
./download-results.py --help
```

!!! tip "No `--profile` necessary on Lightsail for Research"
    Omit `--profile openscpca` from the commands below if you're using Lightsail for Research.

Example script usage below assumes you have created a profile called `openscpca` when [configuring the AWS CLI](../../technical-setup/environment-setup/configure-aws-cli.md).

To download all results for one or more modules, use the `--modules` option as follows:

```sh
# Get results from a single module
./download-results.py --modules module-name --profile openscpca

# Get results from several modules
./download-results.py --modules first-module,second-module --profile openscpca
```

To list all available modules for which you can download results, you can use the `--list-modules` option as follows.
This will list all modules available in, by default, the current release:

```sh
./download-results.py --list-modules --profile openscpca
```

To review what result files would be downloaded without performing the download yet, you can use the `--dryrun` option as follows:

```sh
./download-data.py --modules module-name --dryrun --profile openscpca
```

If you wish to only working with a subset of a module's results, you can use the `--projects` or `--samples` options to download select samples (note: these options are mutually exclusive).
Note that whether this subsetting is possible will depend on whether that module's results have been structured in a manner that reflects this organization. PHRASING HELP PLEASE!

#### Download result data structure

By default, results will download into the `data` directory in your local copy of the `OpenScPCA-analysis` repository, but this can be customized with the `--data-dir` option.

Result downloads will generally follow this structure:

```sh
data
├── {Release}
│       └── results
│          └── {Module name}
│              └── {Result files}
└── current -> {Release}
```

Both result files and ScPCA data files follow the same release schedule.
Therefore, if you have downloaded both data and workflow results, all files will be present in the `data` directory with this overall structure:

```sh
data
├── {Release}
│       ├── {Project ID}
│       │   └── {Sample ID}
│       │       └── {Library files}
│       ├── bulk_metadata.tsv (if applicable)
│       ├── bulk_quant.tsv (if applicable)
│       ├── single_cell_metadata.tsv
|       └── results
|           └── {Module name}
│              └── {Result files}
└── current -> {Release}
```

To list all available releases, you can use the following command from the root of the repository:

```sh
./download-results.py --list-releases
```

The release directory is symlinked to `data/current`.
This is generally the path you should use in your code.
