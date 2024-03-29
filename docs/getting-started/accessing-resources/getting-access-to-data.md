# Getting Access to Data

To participate in the OpenScPCA Project, you will need access to the Single-cell Pediatric Cancer Atlas (ScPCA) data.
Here, we describe several ways to access data.

You can find out more about the contents of files and how they were processed from the ScPCA documentation: <https://scpca.readthedocs.io>.

## Accessing data on the ScPCA Portal

ScPCA data are readily available from the ScPCA Portal that the Data Lab maintains: <https://scpca.alexslemonade.org/>.

You can select the project(s) or sample(s) you are interested in analyzing and download them from the Portal.

## Accessing data on S3

Because we expect that contributors may want to work with many samples or analyze data on remote systems, we also provide access to the ScPCA on AWS S3 via a script for downloading data.

Before you can access data in this manner, the Data Lab team needs to create an AWS account for you; these data are not publicly accessible.
See [Getting Access to Resources](index.md) for more information.

### Using the download data script

!!! note
    These instructions assume you have [configured the AWS CLI](STUB_LINK) or are [using Lightsail for Research](STUB_LINK).

You can list all the options for the download data script by running the following from the root directory of the repository:

```sh
./download-data.py --help
```

The script is designed to download files from whatever release you specify to a folder named for the date of that release and symlink it to `data/current`.



## Accessing simulated data

You can also use the download data script to download the test data for the project by using the `--test-data` option.
The test data are generally smaller, simulated or permuted data files with the same structure as the real project data (learn more at the [`simulate-sce`](https://github.com/AlexsLemonade/OpenScPCA-analysis/tree/main/analyses/simulate-sce) module).

To download test data with all other options set at default, run the following from the root of the repository:

```sh
./download-data.py --test-data
```

This will download test data and direct `data/current` symlink to the test data directory.

To switch back to using ScPCA data, rerun the script with the `--release current` option.