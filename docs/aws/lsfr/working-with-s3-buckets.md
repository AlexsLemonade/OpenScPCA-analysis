# Working with S3 buckets

As an onboarded OpenScPCA contributor, you will be [assigned an AWS S3 results bucket](../../getting-started/accessing-resources/index.md#getting-access-to-aws).
You will sync your analysis results to your bucket [before filing a pull request](../../contributing-to-analyses/creating-pull-requests/index.md).
This is necessary because [analysis module results are not included in Git version control](../../contributing-to-analyses/analysis-modules/index.md#skeleton-analysis-module-contents), so uploading your results to S3 allows the Data Lab to access them during code review.

This page covers some information you'll need for working with your S3 bucket.

!!! note "Still need to join OpenScPCA?"
    Start the process of becoming a contributor by [filing out project intake form](https://share.hsforms.com/1MlLtkGYSQa6j23HY_0fKaw336z0).


## Finding your bucket name

Once you have been accepted as a new OpenScPCA contributor, you will receive an account number specific to you.
Your S3 bucket name is `researcher-{account number}-us-east-2`.

For example, if your account number were `12345`, your bucket name would be `researcher-12345-us-east-2`.

You will need to know your bucket name for a few circumstances:

- You will need to specify your bucket name when [syncing your analysis results with S3](#syncing-your-results-to-s3)
- You will need to tell reviewers your result bucket name when when [filing pull requests (PRs) and filling out the PR template](../../contributing-to-analyses/creating-pull-requests/pull-request-template.md)

### Storing your bucket name

> These instructions are mostly applicable to advanced users, but we're happy to help you set this up if you need assistance!

Optionally, you may wish to save your bucket name as an environment variable in your shell profile file (e.g., `~/.bashrc` for Bash users and `~/.zshrc` for Z shell users).

We have set up some OpenScPCA scripts to look for the environment variable `OPENSCPCA_RESULTS_BUCKET`, so if you define this in your profile, you'll have an easier time using those scripts!

To use this, add this line (with `{account number}` replaced with your account number) to your profile file:

```sh
export OPENSCPCA_RESULTS_BUCKET=researcher-{account number}-us-east-2
```

You can then use `$OPENSCPCA_RESULTS_BUCKET` (or, `${OPENSCPCA_RESULTS_BUCKET}`) to refer to your bucket name when running commands in the terminal or in shell scripts.


## Syncing your results to S3

We have written a script to help you sync your results to your S3 bucket, stored in [`scripts/sync-results.py`](https://github.com/AlexsLemonade/OpenScPCA-analysis/blob/main/scripts/sync-results.py).
For this script to work, you need to be [logged into the AWS account profile you use for contributing to OpenScPCA](../../technical-setup/environment-setup/configure-aws-cli.md#logging-in-to-a-new-session).

This script syncs the contents of a given analysis module's `results` and `plots` directories from your computer to S3.
It does not copy any files from S3 back to your computer.
Your files will be synced to the S3 location `s3://{your bucket name}/{analysis module name}`, where:

- Results will be synced to `s3://{your bucket name}/{analysis module name}/results`
- Plots will be synced to `s3://{your bucket name}/{analysis module name}/plots`


The simplest usage of this script, called from the `OpenScPCA-analysis` repository root folder, is:

```sh
scripts/sync-results.py {analysis module name} \
    --bucket {name of your bucket}
```

- The first argument is the folder name of the analysis module whose results you want to sync
- `--bucket` (or `-b`) is your [bucket name](#finding-your-bucket-name)
  - You can omit this argument if you have saved your bucket name in the environment variable `OPENSCPCA_RESULTS_BUCKET`

**Before you run the script to sync your results,** we recommend you first run the script with the `--dryrun` flag.
This flag will not actually sync your files, but it will instead print out which files will be synced to your S3 bucket, and where they will be placed:

```sh
scripts/sync-results.py {analysis module name} \
    --bucket {name of your researcher bucket} \
    --dryrun
```

If everything looks correct to you when running with `--dryrun`, you can then proceed to run the script without this flag and actually perform the syncing.

### Additional script flags

By default, if there are any result or plot files that exist on S3 but that you have locally deleted, this script will _not delete_ those files from S3.
To override this behavior and delete these files from S3 as well (but only for the specified analysis module), use the `--destructive-sync` flag:

```sh
scripts/sync-results.py {analysis module name} \
    --bucket {name of your researcher bucket} \
    --destructive-sync
```

You can run the following to see all script options:

```sh
scripts/sync-results.py --help
```

#### Advanced options

If you have multiple AWS profiles on your system, it may help to use the `--profile` argument to specify the name of your OpenScPCA AWS profile.

For example, if you [configured your OpenScPCA AWS profile](../../technical-setup/environment-setup/configure-aws-cli.md) to be named `openscpca`, you would use:

```sh
scripts/sync-results.py {analysis module name} \
    --bucket {name of your researcher bucket} \
    --profile openscpca
```
