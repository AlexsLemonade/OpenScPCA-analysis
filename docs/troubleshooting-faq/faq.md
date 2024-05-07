# Frequently asked questions

### Why didn't the sample/project I specified when running the [data download script](../getting-started/accessing-resources/getting-access-to-data.md#using-the-download-data-script) download?

First, we recommend using the `--dryrun` flag when running the script to check which files _would_ be downloaded, to confirm that there is nothing wrong with your internet connection and that you are properly [logged into your AWS profile](../technical-setup/environment-setup/configure-aws-cli.md#logging-in-to-a-new-session).

If you specify a sample or project ID that is not present in the OpenScPCA data release, then nothing will download, but the script will not explicitly fail or warn you about this.

Therefore, you may wish to list all files in the data release S3 bucket to confirm the sample/project files you are trying to download exist.
To do this, you'll need to use the AWS CLI [in a terminal](../software-platforms/general-tools/using-the-terminal.md) and run `aws s3 ls` to list the contents of a given S3 bucket.

All data release files are stored in the S3 bucket `s3://openscpca-data-release/`, which you can recursively list the contents of as follows:

1. First, find the names of all the data release buckets, which are named by release date in [ISO 8601 format](https://en.wikipedia.org/wiki/ISO_8601):
    ```bash
    # List all releases in the overall release bucket
    # If you're _not_ working on Lightsail for Research, specify your AWS profile name with `--profile openscpca`
    aws s3 ls s3://openscpca-data-release/ --profile openscpca
    ```

    !!! tip "Don't forget the trailing slash!"
        The trailing slash at the end of `s3://openscpca-data-release/` is necessary for contents to be listed.
        If you omit the slash and just run `aws s3 ls s3://openscpca-data-release --profile openscpca`, the output will only list the bucket itself, and not its contents..


1. In the output, identify the most recent date.
In the example output below, the most recent release is `2024-05-01/`.
    ```{ .console .no-copy title="Output from listing all releases"}
    PRE 2024-03-11/
    PRE 2024-04-06/
    PRE 2024-05-01/
    ```

1. Then, you can list all projects in the most recent release:
    ```bash
    # List all releases in the 2024-05-01 release bucket, again specifying the AWS profile
    aws s3 ls s3://openscpca-data-release/2024-05-01/ --profile openscpca
    ```

    This has the (abbreviated) output:
    ```{ .console .no-copy title="Output from listing all projects in the 2024-05-01 release"}
    PRE SCPCP000001/
    PRE SCPCP000002/
    PRE SCPCP000003/
    PRE SCPCP000004/
    ```

1. From there, you can continue to find names of nested buckets (all prefixed with `PRE`) and continue listing files.
For example, you can list all files in the `SCPCS000001` sample from the `SCPCP000001` project with:

    ```bash
    # List all SCPCS000001 files in the 2024-05-01 release bucket, again specifying the AWS profile
    aws s3 ls s3://openscpca-data-release/2024-05-01/SCPCP000001/SCPCS000001/ --profile openscpca
    ```

    If you attempt to list a bucket that does not exist, you will not get any output.
    For example, this will not return anything since there is no project ID `SCPCP000099`:
    ```bash
    # List all files in a bucket that does not exist
    aws s3 ls s3://openscpca-data-release/2024-05-01/SCPCP000099/ --profile openscpca
    ```