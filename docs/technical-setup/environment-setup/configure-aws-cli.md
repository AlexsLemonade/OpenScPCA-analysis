# Configure the Amazon Web Services (AWS) Command Line Interface (CLI)

!!! note
    To follow these instructions, we need to create an AWS account for you.
    See [Getting Access to Resources](../../getting-started/accessing-resources/index.md) for more information.

The Amazon Web Services (AWS) Command Line Interface (CLI) allows you to interact with AWS via the command line.
Configuring the AWS CLI is required for using [the download data script](../../getting-started/accessing-resources/getting-access-to-data.md#accessing-data-on-s3) and [uploading your results to S3 for review](STUB_LINK).

These instructions assume you have already installed AWS CLI during [conda set up](setup-conda.md#set-up-conda).
## Configuring the AWS CLI

!!! info
    If you're using [Lightsail for Research](../software-platforms/lsfr/index.md), the configuration and login steps have occurred automatically.
    You only need to follow these steps if you're using another system like your laptop.

To follow these instructions, you'll need [the access portal URL from when you set up your user in IAM Identity Center](../../software-platforms/aws/index.md#joining-iam-identity-center).
Note: SSO below stands for "Single Sign-On" and is what AWS used to call IAM Identity Center.


1. Use the following command in your Terminal application to start configuring:

```sh
aws configure sso
```

2. Fill in the prompts using the following values, replacing `{The access portal URL from you invite email}` including the curly brackets with the access portal URL.

```
SSO session name (Recommended): openscpca-sso
SSO start URL [None]: {The access URL from you invite email}
SSO region [None]: us-east-1
SSO registration scopes [sso:account:access]: sso:account:access
```

3. Step 2 will automatically open a browser window for you to confirm access by hitting the `Confirm and continue` button.
(If a browser window doesn't open, follow the instructions in your Terminal.)
If you're not logged in, you will be prompted to log in.

4. Click the `Allow access` button on the next screen.

5. Return to your Terminal application.
You should see the following:

```
The only AWS account available to you is: {account number}
Using the account ID {account number}
The only role available to you is: ResearcherRestriction
Using the role name "ResearcherRestriction"
```

And fill in the prompts with the following:

```
CLI default client Region [None]: us-east-2
CLI default output format [None]: json
CLI profile name [ResearcherRestriction-{account number}]: openscpca
```

## Logging in to a new session

To use the AWS CLI for this project, you will need to be logged into the `openscpca` profile on your computer.

1. To login, use the following command in your Terminal application:

```sh
aws sso login --profile openscpca
```

2. Step 1 will automatically open a browser window for you to confirm access by hitting the `Confirm and continue` button.
(If a browser window doesn't open, follow the instructions in your Terminal.)

3. Click the `Allow access` button on the next screen.
