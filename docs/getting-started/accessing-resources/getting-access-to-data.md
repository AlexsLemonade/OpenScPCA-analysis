# Getting Access to Data

To participate in the OpenScPCA Project, you will need access to the Single-cell Pediatric Cancer Atlas (ScPCA) data.
Here, we describe several ways to access data.

You can find out more about the contents of files and how they were processed from the ScPCA documentation: <https://scpca.readthedocs.io>.

## Accessing data on the ScPCA Portal

ScPCA data are readily available from the ScPCA Portal that the Data Lab maintains: <https://scpca.alexslemonade.org/>.

You can select the project(s) or sample(s) you are interested in analyzing and download them from the Portal.

We recommend creating a `portal_downloads` subdirectory in the local copy of the `data` directory, which can be accomplished by running the following command from the root of the repository on Linux or Mac OS:

```sh
mkdir -p data/portal_downloads
```

You can then develop your analysis using these paths temporarily.

Before filing a pull request, you should change your paths to reflect the directory established by the download script (`data/current`) described in the next section.

## Accessing data on S3

Because we expect that contributors may want to work with many samples or analyze data on remote systems, we also provide access to the ScPCA on AWS S3 via a script for downloading data.

Before you can access data in this manner, the Data Lab team needs to create an AWS account for you; these data are not publicly accessible.
See [Getting Access to Resources](index.md) for more information.

!!! info Review additional restrictions
    Please review the `DATA_USAGE.md` file included in all downloads for any additional restrictions on data usage (see [Policies](../../policies/index.md)).

### Using the download data script

!!! note
    These instructions assume you have taken the following steps:

    - [Cloned the repo](../../technical-setup/clone-the-repo.md)
    - [Set up your environment](../../technical-setup/environment-setup/index.md),
    - [Configured the AWS CLI and logged in](../../technical-setup/environment-setup/configure-aws-cli.md) OR are [using Lightsail for Research](../software-platforms/lsfr/index.md) (which doesn't require logging in via the AWS CLI)

The `download-data.py` script is designed to download files from whatever release you specify to a folder named for the date of that release and [symlink](https://en.wikipedia.org/wiki/Symbolic_link) it to `data/current`.

We briefly cover some of the most common use cases below, but we encourage you to review all the options available to you.

You can list all the options for the download data script by running the following from the root directory of the repository:

```sh
./download-data.py --help
```

Assuming you created a profile called `openscpca` when [configuring the AWS CLI](../../technical-setup/environment-setup/configure-aws-cli.md) you can run the download script with all default options to download all processed samples in `SingleCellExperiment` format:

```sh
./download-data.py --profile openscpca
```

If you prefer to work with `AnnData` files, you can specify using the `--format` option as follows:

```sh
./download-data.py --format AnnData --profile openscpca
```

To review what samples would be downloaded without performing the download yet, you can use the `--dry-run` option as follows:

```sh
./download-data.py --dry-run --profile openscpca
```

If you're only working with a subset of the data, you can use the `--projects` or `--samples` to download select samples (note: these options are mutually exclusive).

### Download data structure

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
└── current -> {Absolute Path to Repository}/data/{Release}
```

See the ScPCA documentation for more information about individual files: <https://scpca.readthedocs.io>.

Data releases are dated using the following format: `YYYY-MM-DD`.
By default, the data download script will download the most recent release.

The release directory is symlinked to `data/current`.
This is generally the path you should use in your code.

## Accessing simulated data

You can also use the download data script to download the test data for the project by using the `--test-data` option.
The test data are generally smaller, simulated or permuted data files with the same structure as the real project data (learn more at the [`simulate-sce`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/simulate-sce) module).

You do not need an AWS account set up to download the test data.

To download test data with all other options set at default, run the following from the root of the repository:

```sh
./download-data.py --test-data
```

This will download test data and direct `data/current` symlink to the test data directory.

To switch back to using ScPCA data, rerun the script with the `--release current` option.
