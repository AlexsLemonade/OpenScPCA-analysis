# Getting access to data

To participate in the OpenScPCA Project, you will need access to the Single-cell Pediatric Cancer Atlas (ScPCA) data.
Here, we describe several ways to access data.

Broadly speaking, there are three kinds of ScPCA data you might wish to work with:

1. Data from the ScPCA Portal
    - You can find out more about the contents of files and how they were processed from the ScPCA documentation: <https://scpca.readthedocs.io>.
    - OpenScPCA project contributors [with an Amazon Web Services (AWS) account](./index.md#getting-access-to-aws) can access data using the provided `download-data.py` script, [as described below](#using-the-download-data-script).
    - All data is also available [directly from the ScPCA Portal](#accessing-data-from-the-scpca-portal) for users who have not yet fully joined the OpenScPCA project.
2. Results from other OpenScPCA modules
    - OpenScPCA project contributors with an AWS account can obtain results from completed modules using the provided `download-results.py` script, [as described below](#accessing-scpca-module-results).
    - To use results from an in-progress module from the `OpenScPCA-analysis` repository, you may need to run the module yourself to generate its result files.
3. Test datasets and associated results
    - We provide reduced-size test files with simulated and/or permuted versions of both ScPCA Portal data and results from completed modules.
    - These data are used for automated testing, but you can also use them while developing your module if smaller files helps make development more efficient.
    - You do not need an AWS account to download the test files, so you can use this data before you are granted full access to ScPCA data.


## Accessing data as an OpenScPCA contributor

Because we expect that contributors may want to work with many samples or analyze data on remote systems, contributors can access to ScPCA data and results from AWS S3.
We provide a download script to make accessing these data more convenient.

Before you can access data in this manner, the Data Lab team needs to create an AWS account for you, as these data are not publicly accessible.
See our documentation [getting access to AWS](index.md#getting-access-to-aws) for more information.

!!! info Review additional restrictions
    Please review the `DATA_USAGE.md` file included in all downloads for any additional restrictions on data usage (see [Policies](../../policies/index.md)).

### Using the download data script

TODO DURING REVIEW: Do we still need this admonition note?
!!! note
    These instructions assume you have taken the following steps:

    - [Cloned the repo](../../technical-setup/clone-the-repo.md)
    - [Set up conda](../../technical-setup/environment-setup/setup-conda.md) (note that conda is pre-installed, but not yet set up, on [Lightsail for Research](../../aws/lsfr/starting-development-on-lsfr.md#create-and-activate-a-conda-environment) instances)
    - [Configured the AWS CLI and logged in](../../technical-setup/environment-setup/configure-aws-cli.md) OR are [using Lightsail for Research](../../aws/index.md#lightsail-for-research-virtual-computing-with-aws) (which doesn't require logging in via the AWS CLI)

The [`download-data.py` script](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/download-data.py) is designed to download files from whatever release you specify to a directory in `data` named for the date of that release and [symlink](https://en.wikipedia.org/wiki/Symbolic_link) it to `data/current`.

We briefly cover some of this script's most common use cases below, but we encourage you to review all the options available to you.
All examples below assume you are running the script from the root of the repository.

!!! tip "Log into AWS CLI before running the script"
    Before running this script locally (i.e., not using [Lightsail for Research](../../aws/index.md#lightsail-for-research-virtual-computing-with-aws)), you will need to be [logged into your AWS CLI profile](../../technical-setup/environment-setup/configure-aws-cli.md#logging-in-to-a-new-session).
    To do this, run the following commands in terminal before running the `download-data.py` script, and follow instructions to log in.

    ```sh
    # replace `openscpca` with your AWS CLI profile name if it differs
    export AWS_PROFILE=openscpca
    aws sso login
    ```

    The command `export AWS_PROFILE=openscpca` will define your AWS profile name as `openscpca` for the duration of your terminal session.
    Defining this environment variable is helpful because the `download-data.py` script also needs your profile name to download data from S3.

    Note that you can also [add this profile definition to your shell profile file](../../technical-setup/environment-setup/configure-aws-cli.md#storing-your-aws-profile-name) so that it will always be defined in any terminal session.

- You can list all available arguments and their usage for the download data script with the following:

    ```sh
    ./download-data.py --help
    ```

- You can run the download script with all default options to download all processed samples from the most recent release in `SingleCellExperiment` format:

    !!! warning
        **The full data release download is quite large: over 35 GB for `SingleCellExperiment` format, and over 100 GB for `AnnData` format.
        Please download the full data release with caution!**

    ```sh
    ./download-data.py
    ```

- If you prefer to work with `AnnData` files, use the `--format` option as follows:

    ```sh
    ./download-data.py --format AnnData
    ```

- To download a specific release other than the most recent release, you can use the `--release` option.
Again, download full data releases with caution, as they are quite large!
    - You can use the `--list-releases` flag to see all available releases, which are named in `YYYY-MM-DD` format based on their release date.

    ```sh
    # First, see which releases are available for download
    ./download-data.py --list-releases

    # Download SingleCellExperiment format files for a specific release
    # For example, this downloads the 2024-05-01 release
    ./download-data.py --release 2024-05-01

    # Download AnnData format files for a specific release, for example
    ./download-data.py --release 2024-05-01 --format AnnData
    ```

- If you're only working with a subset of the data, you can use the `--projects` or `--samples` option to download select projects or samples, respectively (note that these options are mutually exclusive):

    ```sh
    # use the --projects option
    ./download-data.py --projects SCPCPXXXXXX

    # use the --samples option
    ./download-data.py --samples SCPCSXXXXXX,SCPCSXXXXXY
    ```

- To review the files would be downloaded without performing the download yet, you can use the `--dryrun` option as follows.
    - We strongly recommend using the `--dryrun` flag the first time you run the script for a new data download to ensure the downloaded files are as expected.

    ```sh
    ./download-data.py --dryrun

    # show what would be downloaded when the --projects option is used, for example
    ./download-data.py --projects SCPCPXXXXXX --dryrun
    ```

The `download-data.py` script will update the `data/current` symlink to point to the current release directory when run without a `--release` option or to the test data directory when run with `--test-data`.
If you want to instead point to a specific release, you can use the [`--update-symlink` flag described below](#updating-the-current-symlink) to update which release directory the `data/current` symlink directs to.


### Downloaded data file structure

Downloads via the download script will generally have the following structure.

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

!!! tip "The `data/current` directory"
    The `data/current` directory is a [symlink](https://en.wikipedia.org/wiki/Symbolic_link) to the most recent data release.
    You should therefore use the `data/current` path in your analysis code when reading in data, rather than a specific release directory.


See the ScPCA documentation for [more information about individual files](https://scpca.readthedocs.io/en/latest/sce_file_contents.html).

Data releases are dated using the following format: `YYYY-MM-DD`.
By default, the data download script will download the most recent release.

To list all available releases, you can use the following command from the root of the repository:

```sh
./download-data.py --list-releases
```

## Accessing ScPCA module results

As an analysis modules matures, the Data Lab will port the analysis to a separate repository called [`OpenScPCA-nf`](https://github.com/AlexsLemonade/OpenScPCA-nf).
This repository holds a [Nextflow workflow](https://www.nextflow.io/) to reproducibly run and generate results from completed analysis modules.

If you wish to use output from an `OpenScPCA-nf` analysis module in your module, you will need to obtain the results generated by the workflow.

We have provided a [`download-results.py` script](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/download-results.py) that you can use to download results from a module that has been run in `OpenScPCA-nf`.

We briefly cover some of the most common use cases below, but we encourage you to review all the options available to you.
All examples below assume you are running the script from the root of the repository.

!!! tip "Log into AWS CLI before running the script"
    Before running this script locally (i.e., not from an [ALSF-provided virtual computer](../../aws/index.md#lightsail-for-research-virtual-computing-with-aws)), you will need to be [logged into your AWS CLI profile](../../technical-setup/environment-setup/configure-aws-cli.md#logging-in-to-a-new-session).
    To do this, run the following commands in terminal before running the `download-results.py` script, and follow instructions to log in.

    ```sh
    # replace `openscpca` with your AWS CLI profile name if it differs
    export AWS_PROFILE=openscpca
    aws sso login
    ```

    The command `export AWS_PROFILE=openscpca` will define your AWS profile name as `openscpca` for the duration of your terminal session.
    Defining this variable is helpful because the `download-data.py` script also needs your profile name to download data from S3.

    Note that you can also [add this profile definition to your shell profile file](../../technical-setup/environment-setup/configure-aws-cli.md#storing-your-aws-profile-name) so that it will always be defined in any terminal session.

- You can list all available arguments and their usage for the results download script with the following:

    ```sh
    ./download-results.py --help
    ```

- To download results for one or more modules, use the `--modules` option with the module directory name(s) as a comma separated list:

    ```sh
    # Get results from a single module
    ./download-results.py --modules module-name

    # Get results from several modules
    ./download-results.py --modules first-module,second-module
    ```

- To list all available modules for which you can download results, you can use the `--list-modules` flag as follows.
This will list all modules available in, by default, the current release:

    ```sh
    ./download-results.py --list-modules
    ```

- To see the result files that would be downloaded without performing the download yet, you can use the `--dryrun` flag as follows:
    - We strongly recommend using the `--dryrun` flag the first time you run the script for a new data download to ensure the downloaded files are as expected.

    ```sh
    ./download-results.py --modules module-name --dryrun
    ```

- While the structure of results files varies by module, some modules' results will be organized by project and/or sample.
For those modules, you can use either the `--projects` or `--samples` flag (noting that these options are mutually exclusive) to download results specific to one or more projects or samples, respectively:

    ```sh
    # use the --projects option
    ./download-results.py --modules module-name --dryrun --projects SCPCPXXXXXX

    # use the --samples option
    ./download-results.py --modules module-name --dryrun --samples SCPCSXXXXXX,SCPCSXXXXXY
    ```


### Downloaded results file structure

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

!!! tip "The `data/current` directory"
    The `data/current` directory is a [symlink](https://en.wikipedia.org/wiki/Symbolic_link) to the most recent release directory available locally.
    You should therefore use the `data/current` path in your analysis code when reading in results, rather than a specific release directory.


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

## Accessing test data

The test data are simulated or permuted data files with the same structure as the real project data, and are generally much smaller than the original data (learn more at the [`simulate-sce`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/simulate-sce) module).
These files are also used during automated testing of analysis modules. <!-- STUB_LINK for module GHAs -->

You can use the [download data script](#using-the-download-data-script) to download the test data for the project by using the `--test-data` flag.
Similarly, you use the [download result script](#accessing-scpca-module-results) to download ScPCA results from running completed modules on the test data by using the `--test-data` flag.

You do not need an AWS account set up to download the test data or results.

- To download test data with all other options set at default, run the following from the root of the repository:

    ```sh
    ./download-data.py --test-data
    ```

- To download a given module's results as generated with the test data, run the following from the root of the repository:

    ```sh
    ./download-results.py --modules module-name --test-data
    ```


Running either of the above commands will update the `data/current` symlink to point to the `data/test` directory.
This means any real ScPCA data or results you had previously downloaded will no longer be in the `data/current` path.
You can use the `--update-symlink` flag, as described below, to update which release directory the `data/current` symlink directs to.

## Updating the `current` symlink

As described above, the `data/current` directory is a [symlink](https://en.wikipedia.org/wiki/Symbolic_link) to the most recent release directory.
You can update this symlink to point to a different data release directory (or to the test data directory) in one of two ways:

1. When you download a new data or results release, `data/current` will automatically be updated to direct to the newly-downloaded release directory.
    - Both the `download-data.py` and `download-results.py` scripts will print a message telling you where `current` directs to.
1. You can use the `download-data.py` script with the `--update-symlink` flag to redirect the `data/current` symlink to a release of your choice.
Note that no data will be downloaded if you use this flag.
Below are some common use cases for this flag:

    ```sh
    # By default, this flag will update the current symlink to direct to the
    #  most recent release that is present on your computer
    ./download-data.py --update-symlink

    # Use the --release option to specify which release to direct the symlink to
    ./download-data.py --update-symlink --release 2025-05-01

    # Use the --test-data flag to direct the symlink to the test data
    ./download-data.py --update-symlink --test-data
    ```


## Accessing data from the ScPCA Portal

If you have not formally joined the OpenScPCA project, you will need to obtain data directly from ScPCA Portal that the Data Lab maintains: <https://scpca.alexslemonade.org/>.

You can select the project(s) or sample(s) you are interested in analyzing and download them from the Portal.

We recommend creating a `portal-downloads` subdirectory in the local copy of the `data` directory, which can be accomplished by running the following command from the root of the repository on Linux or macOS:

```sh
mkdir -p data/portal-downloads
# set up a symlink
ln -s data/portal-downloads data/current
```

You can then develop your analysis using these paths.

However, before filing a pull request, you should change your paths to use `data/current` directory [established by the data download script](#downloaded-data-file-structure).