# Getting access to compute

## Lightsail for Research

We use [Lightsail for Research](https://aws.amazon.com/lightsail/research/) to make virtual computers available to OpenScPCA contributors.

You must meet the criteria for AWS account creation to gain access to Lightsail for Research.
See [Getting access to AWS](index.md#getting-access-to-aws) for more information.

See our docs about [Lightsail for Research](../../aws/index.md#lightsail-for-research-virtual-computing-with-aws) for more information on using this service.

### Monthly budget

All contributors have a 200 USD/month budget for compute usage that they must monitor.

Lightsail for Research virtual computers and volumes will report their monthly cost.
See [the AWS documentation on monitoring cost estimates](https://docs.aws.amazon.com/lightsail-for-research/latest/ug/monitor-cost-usage-estimates.html) for more information.

[Your results S3 bucket](../../aws/working-with-s3-buckets.md) also counts towards your budget.
You can not access the cost estimates for S3 directly, so please use $0.023/GB as your as a rule of thumb for [S3 storage pricing in `us-east-2`](https://aws.amazon.com/s3/pricing/).

When you have exhausted your budget, we may prevent you from accessing services for the rest of the month.
We will notify you via the email you used to sign up that you are approaching your budget before restricting access.

Requests for an increased budget can be included on Discussions or Issues in sections that ask about computing resources.

Requests for a budget increase will be evaluated on a case-by-case basis.

### GPU instance access

Accessing GPU instances in Lightsail for Research in your account requires us to put in a request with Amazon Web Services.

Requests for access to GPU instances can be included on Discussions or Issues in sections that ask about computing resources.

Requests for access to GPU instances will be evaluated on a case-by-case basis.
