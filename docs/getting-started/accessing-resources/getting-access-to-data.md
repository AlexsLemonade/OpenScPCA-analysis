# Getting access to data

To participate in the OpenScPCA Project, you will need access to the Single-cell Pediatric Cancer Atlas (ScPCA) data.
Here, we describe several ways to access data.

Broadly speaking, there are three kinds of ScPCA data you might wish to work with:


- Data from the ScPCA Portal
    - You can find out more about the contents of files and how they were processed from the ScPCA documentation: <https://scpca.readthedocs.io>
- Results from other OpenScPCA modules
    - You can obtain results from modules that have been added to the `OpenScPCA-nf` workflow by [following the instructions below](#accessing-scpca-module-results).
    - To use results from an in-progress module from the `OpenScPCA-analysis` repository, you may need to run the module yourself to generate its result files.
- Test datasets
    - We provide reduced-size test files with simulated and/or permuted versions of both ScPCA Portal data and results from completed modules.
    - These data are used for automated testing, but you can also use them while developing your module if smaller files helps make development more efficient.
    - You do not need an AWS account to download the test files, so you can use this data before you are granted full access to ScPCA data.

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

Because we expect that contributors may want to work with many samples or analyze data on remote systems, access to ScPCA data on AWS S3 is also available to contributors.
We also provide a download script to make accessing these data more convenient.

Before you can access data in this manner, the Data Lab team needs to create an AWS account for you; these data are not publicly accessible.
See our documentation [getting access to AWS](index.md#getting-access-to-aws) for more information.

!!! info Review additional restrictions
    Please review the `DATA_USAGE.md` file included in all downloads for any additional restrictions on data usage (see [Policies](../../policies/index.md)).

#### Using the download data script

!!! note
    These instructions assume you have taken the following steps:

    - [Cloned the repo](../../technical-setup/clone-the-repo.md)
    - [Set up conda](../../technical-setup/environment-setup/setup-conda.md) (note that conda is pre-installed, but not yet set up, on [Lightsail for Research](../../aws/lsfr/starting-development-on-lsfr.md#create-and-activate-a-conda-environment) instances)
    - [Configured the AWS CLI and logged in](../../technical-setup/environment-setup/configure-aws-cli.md) OR are [using Lightsail for Research](../../aws/index.md#lightsail-for-research-virtual-computing-with-aws) (which doesn't require logging in via the AWS CLI)

The [`download-data.py` script](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/download-data.py) is designed to download files from whatever release you specify to a folder named for the date of that release and [symlink](https://en.wikipedia.org/wiki/Symbolic_link) it to `data/current`.

!!! tip "Log into AWS CLI before running the script"
    Before running this script locally, you will need to be [logged into your AWS CLI profile](../../technical-setup/environment-setup/configure-aws-cli.md#logging-in-to-a-new-session).
    To do this, run the following commands in terminal before running the `data-download.py` script, and follow instructions to log in.

    ```sh
    # replace `openscpca` with your AWS CLI profile name if it differs
    export AWS_PROFILE=openscpca
    aws sso login
    ```

    The command `export AWS_PROFILE=openscpca` will define your AWS profile name as `openscpca` for the duration of your terminal session.
    Defining this environment variable is helpful because the `download-data.py` script also needs your profile name to download data from S3.

    Note that you can also [add this profile definition to your shell profile file](../../technical-setup/environment-setup/configure-aws-cli.md#storing-your-aws-profile-name) so that it will always be defined in any terminal session.

We briefly cover some of the most common use cases below, but we encourage you to review all the options available to you.

You can list all the options for the download data script by running the following from the root directory of the repository:

```sh
./download-data.py --help
```

You can run the download script with all default options to download all processed samples from the most recent release in `SingleCellExperiment` format:

!!! warning
    **The full data release download is quite large: over 35 GB for `SingleCellExperiment` format, and over 100 GB for `AnnData` format.
    Please download the full data release with caution!**

```sh
./download-data.py
```

If you prefer to work with `AnnData` files, use the `--format` option as follows:

```sh
./download-data.py --format AnnData
```

To review what samples would be downloaded without performing the download yet, you can use the `--dryrun` option as follows:

```sh
./download-data.py --dryrun
```

If you're only working with a subset of the data, you can use the `--projects` or `--samples` to download select samples (note: these options are mutually exclusive):

```sh
# use the --projects option
./download-data.py --projects SCPCPXXXXXX

# use the --samples option
./download-data.py --samples SCPCSXXXXXX,SCPCSXXXXXY
```

#### Download data structure

Downloads via the download script will generally have the following structure:

```sh
data
├── {Release}
│       └── {Project ID}
│           ├── {Sample ID}
│           │     └── {Library files}
│           ├── bulk_metadata.tsv (if applicable)
│           ├── bulk_quant.tsv (if applicable)
│           └── single_cell_metadata.tsv
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


## Accessing ScPCA module results

As an analysis modules matures, the Data Lab will port the analysis to a separate repository called [`OpenScPCA-nf`](https://github.com/AlexsLemonade/OpenScPCA-nf).
This repository [holds a Nextflow workflow](../../ensuring-repro/openscpca-nf/index.md) to reproducibly run and generate results from completed analysis modules.

If you wish to use output from an `OpenScPCA-nf` analysis module in your module, you will need to obtain the results generated by the workflow.

We have provided a [`download-results.py` script](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/download-results.py) that you can use to download results from a module that has been run in `OpenScPCA-nf`.

We briefly cover some of the most common use cases below, but we encourage you to review all the options available to you.

You can list all the options for `download-results.py` by running the following from the root directory of the repository:

```sh
./download-results.py --help
```

!!! tip "Log into AWS CLI before running the script"
    Before running this script locally, you will need to be [logged into your AWS CLI profile](../../technical-setup/environment-setup/configure-aws-cli.md#logging-in-to-a-new-session).
    To do this, run the following commands in terminal before running the `results-download.py` script, and follow instructions to log in.

    ```sh
    # replace `openscpca` with your AWS CLI profile name if it differs
    export AWS_PROFILE=openscpca
    aws sso login
    ```

    The command `export AWS_PROFILE=openscpca` will define your AWS profile name as `openscpca` for the duration of your terminal session.
    Defining this variable is helpful because the `download-data.py` script also needs your profile name to download data from S3.

    Note that you can also [add this profile definition to your shell profile file](../../technical-setup/environment-setup/configure-aws-cli.md#storing-your-aws-profile-name) so that it will always be defined in any terminal session.

To download results for one or more modules, use the `--modules` option as follows:

```sh
# Get results from a single module
./download-results.py --modules module-name

# Get results from several modules
./download-results.py --modules first-module,second-module
```

To list all available modules for which you can download results, you can use the `--list-modules` option as follows.
This will list all modules available in, by default, the current release:

```sh
./download-results.py --list-modules
```

To see the result files that would be downloaded without performing the download yet, you can use the `--dryrun` option as follows:

```sh
./download-results.py --modules module-name --dryrun
```

While the structure of results files varies by module, some modules' results will be organized by project and/or sample.
For those modules, you can use either the `--projects` or `--samples` flag (noting that these options are mutually exclusive) to download results specific to one or more projects or samples:

```sh
# use the --projects option
./download-results.py --modules module-name --dryrun --projects SCPCPXXXXXX

# use the --samples option
./download-results.py --modules module-name --dryrun --samples SCPCSXXXXXX,SCPCSXXXXXY
```


#### Download result data structure

Results will download into the `data` directory in your local copy of the `OpenScPCA-analysis` repository.
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
│       │   │     └── {Library files}
│       │   ├── bulk_metadata.tsv (if applicable)
│       │   ├── bulk_quant.tsv (if applicable)
│       │   └── single_cell_metadata.tsv
|       └── results
|           └── {Module name}
│              └── {Result files}
└── current -> {Release}
```

To list all available releases, you can use the following command from the root of the repository:

```sh
./download-results.py --list-releases
```

The `data/current` directory is symlinked to the current release directory
.
This is generally the path you should use in your code.


## Accessing test data

The test data are simulated or permuted data files with the same structure as the real project data, and are generally much smaller than the original data (learn more at the [`simulate-sce`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/simulate-sce) module).
These files are also used during [automated testing of analysis modules](../../ensuring-repro/workflows/run-module-gha.md).

You can use the [download data script](#using-the-download-data-script) to download the test data for the project by using the `--test-data` option.
Similarly, you use the [download result script](#accessing-scpca-module-results) to download ScPCA results from running completed modules on the test data by using the `--test-data` option.

You do not need an AWS account set up to download the test data or results.

To download test data with all other options set at default, run the following from the root of the repository:

```sh
./download-data.py --test-data
```

To download a given module's results as generated with the test data, run the following from the root of the repository:

```sh
./download-results.py --modules module-name --test-data
```

<!--

TODO: This section needs to be reworked to avoid unexpected downloads. We should recommend --use-release instead.

These commands will download either test data or results, respectively, and will update the `data/current` symlink to instead point to the `data/test` directory.
This means any real ScPCA data or results you had previously downloaded will no longer be in the `data/current` path.
To switch this symlink path back to ScPCA data or results, rerun either script with the `--release current` option:

```sh
./download-data.py --release current

./download-results.py --release current
```
If you had already downloaded the most recent data or results, this will not repeat downloading the files you already have.

-->