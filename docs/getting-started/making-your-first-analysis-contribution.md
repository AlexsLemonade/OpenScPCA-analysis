# Making your first analysis contribution

<figure markdown="span">
      ![Path from joining to making a contribution](../img/first-analysis-contribution.png){width="700"}
</figure>

This is the path that you can expect to follow when making your first contribution to an analysis.

* [Planning](#planning)
* [Implementation](#implementation)
* [Review](#review)

## Planning

### Make an analysis proposal

[Start with GitHub Discussions](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/111) to propose a new analysis or a modification to an existing analysis.

You will file a Discussion for:

* A new analysis in the [`Propose a new analysis` category](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/categories/propose-a-new-analysis)
* An update to an existing analysis in the [`Modify an existing analysis` category](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/categories/modify-an-existing-analysis)

A Data Lab member will review your proposal and respond to your original post to discuss further.
Together you will continue to develop your analysis plans and make decisions about next steps.

### File an issue

* After you have proposed an analysis on Discussions and a plan has been developed, then you can file [an issue](../communications-tools/github-issues/index.md) in the `AlexsLemonade/OpenScPCA-analysis` repository.
* Refer to the documentation on [how to write an issue](../communications-tools/github-issues/writing-issues.md).
* We use [issue templates](../communications-tools/github-issues/issue-templates.md) to prompt issue authors to [include helpful information](../communications-tools/github-issues/what-makes-a-good-issue.md) that helps us accomplish the work tracked in that Issue.

### Pick an existing issue

* There may be existing issues in the `AlexsLemonade/OpenScPCA-analysis` repository for planned analyses that have not been started yet.
* If you would like to take on a planned analysis, you can comment on the issue to note your interest, ask any clarifying questions, and/or propose a solution for the Data Lab team to review and discuss further.

### Getting access to resources

* Once you've met [the criteria for Amazon Web Services (AWS) account creation](accessing-resources/index.md#getting-access-toaws), which includes discussing an analysis, we will create an AWS account for you.
* You will need to accept the invitation to join [as described in our AWS documentation](../software-platforms/aws/index.md#joining-iam-identity-center) to complete one of the steps in environment setup mentioned below.

## Implementation

### Create local setup

Before starting your analysis, you will need to [set up your local environment](../technical-setup/index.md).

* Download and set up a [Git client](../technical-setup/install-a-git-client.md)
* [Fork the `AlexsLemonade/OpenScPCA-analysis` repository](../technical-setup/fork-the-repo.md)
* [Clone your fork](../technical-setup/clone-the-repo.md) to your computer
* Set up [additional dependencies](../technical-setup/environment-setup/index.md) on your computer that you'll need to contribute to OpenScPCA

### Scope work

* When working on an analysis, it's important to take some time to plan your work.
* Refer to this documentation on [scoping your first pull request](../contributing-to-analyses/creating-pull-requests/scoping-pull-requests.md).

### Implement your analysis

* Now that you have set up your local environment, discussed your proposed analysis, and filed an Issue, you're ready to create an [analysis module](../contributing-to-analyses/analysis-modules/index.md).

* Refer to the documentation on [creating an analysis module](../contributing-to-analyses/analysis-modules/creating-a-module.md) to make your first analysis folder.

* [Document each step of your analysis](../contributing-to-analyses/analysis-modules/documenting-analysis.md).

* Once you have added or modified code for your analyses that you are happy with, you will need to [commit your changes](../contributing-to-analyses/working-with-git/making-commits.md) to your [feature branch](../contributing-to-analyses/working-with-git/working-with-branches.md).


## Review

### File a pull request

* Once you have committed changes to your feature branch and [pushed to origin](../contributing-to-analyses/working-with-git/push-to-origin.md), you are ready to file a [pull request](../contributing-to-analyses/creating-pull-requests/index.md).


* Similar to issue templates, we use [pull request templates](../contributing-to-analyses/creating-pull-requests/pull-request-template.md) to prompt you for information that will help us review your code.

### Undergo analytical code review

* [Pull requests must be reviewed and approved](../contributing-to-analyses/creating-pull-requests/index.md#the-pull-request-review-process) by a Data Lab member before they can be merged into the `AlexsLemonade/OpenScPCA-analysis` repository.
* All pull requests will undergo [analytical code review](../contributing-to-analyses/pr-review-and-merge/index.md).

### Merge

* Once approved, a Data Lab staff member will merge your feature branch into the main branch of `AlexsLemonade/OpenScPCA-analysis`. <!-- STUB_LINK: add link to merge docs -->
