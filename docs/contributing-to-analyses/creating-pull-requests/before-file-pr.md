# Before you file a pull request

Before you are able to file a pull request with proposed changes, you will need to complete the following steps.

- Merge any changes from `AlexsLemonade/OpenScPCA-analysis:main` into your [feature branch](../working-with-git/working-with-branches.md).
- Add any results files needed for review to S3.
- Determine software and compute requirements.

## Merge changes from upstream repository

In the time that you are working on your analysis on your [feature branch](../working-with-git/working-with-branches.md), it's likely that others have also been working on changes that have been merged into the `main` branch of `AlexsLemonade/OpenScPCA-analysis`.
You will need to make sure any changes made to the `main` branch are incorporated into your feature branch before filing a PR.

You can do this in GitKraken by finding the remote copy of the `main` branch in the branch graph (indicated by presence of the ALSF lemon logo).
Right-click on the `main` branch and chooose `merge AlexsLemonade/main into username/name-of-feature-branch`.

<figure markdown="span">
    ![Merge main](../../img/before-file-1.png){width="600"}
</figure>

Aftr merging, be sure to [push to origin](../working-with-git/push-to-origin.md), ensuring the changes are reflected in the remote copy of your feature branch.

!!! tip

    If GitKraken notifies you of merge conflicts, see our guide on [resolving merge conflicts](./resolve-merge-conflicts.md).

## Add results files to S3

As a contributor to OpenScPCA, you will be provided with an [AWS account](../../software-platforms/aws/joining-aws.md) and access to an [AWS S3 results bucket](../../software-platforms/aws/working-with-s3-buckets.md).

Results files from running your analysis should not be [under version control and should not be included in your PR](../analysis-modules/index.md#skeleton-analysis-module-contents).
All results files should be saved in your S3 bucket.
Adding your results to S3 allows the Data Lab to access them during code review.

Follow [these instructions to sync your results to S3](../../software-platforms/aws/working-with-s3-buckets.md#syncing-your-results-to-s3) prior to filing your PR.

## Determine software and compute requirements

When filling out the [PR template](./pull-request-template.md) you will be asked to provide information on the [computational resources and software requirements](../determining-requirements/index.md) needed for running your proposed analysis.

See the [documentation on software requirements](../determining-requirements/determining-software-requirements.md) to identify how to track and determine any software requirements.
Please include any files used to track software requirements, such as `renv.lock` or `environment.yml`, in your PR.

See the [documentation on computational resources](../determining-requirements/determining-compute-requirements.md) to identify the computational resources needed for your analysis.
