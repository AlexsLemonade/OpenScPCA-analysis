# Setting up AWS: Joining IAM Identity Center

When we create an AWS account for you, we will also create a user for you with the username `researcher-{your github handle}`.

You will receive an email from `no-reply@login.awsapps.com` with the subject "Invitation to join AWS IAM Identity Center (successor to AWS Single Sign-On)."

You must accept this invitation to gain access to your AWS account and to download project data from S3.

To accept the invitation, please take the following steps.

1. Note the following information in the invitation email: **Your AWS access portal URL**.
   You will need this value to [configure AWS CLI locally](../../technical-setup/environment-setup/configure-aws-cli.md) and to [log into the console to use Lightsail for Research](STUB_LINK).
   It is also helpful to bookmark this URL so you can easily return to it.

2. Click the `Accept Invitation` button in the email.

3. Enter a strong password that is unique to this service when prompted in the browser.

4. Set up a multifactor authentication device, such as an authenticator app or security key, when prompted in the browser.
