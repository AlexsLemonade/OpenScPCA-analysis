# Cloud storage and compute

We use Amazon Web Services (AWS) to grant contributors access to data and Linux virtual computers.

!!! note
    See [Getting Access to Resources](../getting-started/accessing-resources/index.md) for more information about getting AWS access.

    Once we have created an AWS account for you, you can proceed to [set up your account](./joining-aws.md).

    Note that your use of AWS services described in this section counts towards your [monthly budget](../getting-started/accessing-resources/getting-access-to-compute.md#monthly-budget).

## S3: Data storage with AWS

We use [AWS S3](https://aws.amazon.com/s3/) to store project data in folders referred to in S3 as "buckets."
We also use S3 as central location for you to share your result files with us before [filing a pull request](../contributing-to-analyses/creating-pull-requests/index.md).
All contributors will have a designated [researcher bucket](working-with-s3-buckets.md) where results can be synced and shared as part of [code review](../contributing-to-analyses/pr-review-and-merge/index.md).

Please refer to these pages about working with S3:

- [Configuring and logging into AWS from the command line](../technical-setup/environment-setup/configure-aws-cli.md)
- [Accessing project data from S3](../getting-started/accessing-resources/getting-access-to-data.md#accessing-data-from-s3)
- [Using your researcher bucket](working-with-s3-buckets.md)

## Lightsail for Research: Virtual computing with AWS

As an OpenScPCA contributor, you are eligible for a monthly allocation on [Lightsail for Research (LSfR)](https://aws.amazon.com/lightsail/research/), an Amazon product that allows you to work on a virtual Ubuntu (Linux) Desktop computer in your browser.

Please refer to these pages to learn more about setting up and using LSfR:

- [Creating virtual computers](./lsfr/creating-vcs.md)
- [Accessing virtual computers](./lsfr/accessing-vcs.md)
- [Working with storage volumes](./lsfr/working-with-volumes.md)
- [Developing on virtual computers](./lsfr/starting-development-on-lsfr.md)
    - This page also provides instructions for setting up virtual computers after they've been created.
- [Working with snapshots](./lsfr/working-with-snapshots.md)
