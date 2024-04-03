# Working with S3 buckets

As an onboarded OpenScPCA contributor, you will be [assigned an AWS S3 bucket](../../getting-started/accessing-resources/index.md#getting-access-to-aws).
You will use this bucket to sync your analysis results.
This is necessary because [analysis module results are not included in Git version control](../../contributing-to-analyses/analysis-modules/index.md#skeleton-analysis-module-contents), so uploading your results to S3 allows the Data Lab to access them during code review.

This page covers some information you'll need for working with your S3 bucket.

!!! note "Still need to join OpenScPCA?"
    Start the process of becoming a contributor by [filing out project intake form](https://share.hsforms.com/1MlLtkGYSQa6j23HY_0fKaw336z0).


## Finding your bucket name

Once you have been accepted as a new OpenScPCA contributor, you will receive an account number specific to you.
Your S3 bucket name is simply `researcher-{account number}-us-east-2`.

For example, if your account number were `12345`, your bucket name would be `researcher-12345-us-east-2`.

You will need to know your bucket name for a few circumstances:

- You will need to specify your bucket name when [syncing your analysis results with S3](#syncing-your-results-to-s3)
- You will need to tell reviewers your result bucket name when when [filing pull requests (PRs) and filling out the PR template](../../contributing-to-analyses/creating-pull-requests/pull-request-template.md)

### Storing your bucket name

> These instructions are mostly applicable to advanced users, but we're happy to help you set this up if you need assistance!

Optionally, you may wish to save your bucket name as an environment variable in your shell profile file (e.g., `~/.bashrc` for Bash users and `~/.zshrc` for Z shell users).

We have set up some OpenScPCA scripts to look for the environment variable `OPENSCPCA_RESULTS_BUCKET`, so if you define this in your profile, you'll have an easier time using those scripts!

For this, add this line (with `{account number}` replaced with your account number) to your profile file:
```sh
export OPENSCPCA_RESULTS_BUCKET=researcher-{account number}-us-east-2
```

You can then use `$OPENSCPCA_RESULTS_BUCKET` (or, `${OPENSCPCA_RESULTS_BUCKET}`) to refer to your bucket name when running commands in terminal.


## Syncing your results to S3

_Coming up in the next PR._
