# Working with branches

In OpenScPCA we use branches to make any changes to the main code base.
This allows contributors to work on their changes in an isolated environment.
The contributor can then request for changes from that branch to be incorporated via a [pull request](STUB-LINK for pull requests).

## What is a branch?

A branch in GitHub refers to an isolated copy of the repository contents.
For example, the `main` branch of the `OpenScPCA-analysis` repository contains the main copy of all items in the repository.

!!! note "Learn more about branches"
    For more details on branches, see:

    - [GitHub's documentation describing branches](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-branches)
    - [Tutorial on creating branches with GitKraken](https://www.gitkraken.com/learn/git/problems/create-git-branch)

## What is a feature branch?

Instead of directly making code changes in the `main` branch, we write code in a _feature branch_ created off of the `main` branch.
A feature branch contains a copy of the code in an isolated space.

By using a feature branch, you can make code changes without affecting the code in the `main` branch or interfering with other contributors as they make changes.

### When should you use a feature branch?

Every time you make a new addition or modify existing code, you should create a _new_ feature branch.

Every pull request should contain changes made in a _single_ feature branch.
This means you will create a new feature branch for every new pull request you plan to file.

## Overview of working with branches in OpenScPCA

- When you are working on your own analyses or proposing a change to the main code base, create a feature branch from the `main` branch of the `AlexsLemonade/OpenScPCA-analysis` repository.
    - Haven't forked yet? Please see the documentation on [how to create your fork](../../technical-setup/fork-the-repo.md).
    - Be sure that you have added the [`AlexsLemonade/OpenScPCA-analysis` as the upstream repository](../../technical-setup/clone-the-repo.md#add-the-project-repository-as-a-remote-repository).
- [Commit](STUB_LINK to making commits) any changes that you make to that feature branch.
- [Push](STUB_LINK for pushing to origin) your commits to your origin repository, i.e. your fork on GitHub.com
- File a pull request to the `main` branch of the upstream `AlexsLemonade/OpenScPCA-analysis` repository.
  - A [pull request](STUB_LINK for section on pull request) is a request to add code in your feature branch to the main code base.

## Creating a feature branch in GitKraken

Follow these steps to create your feature branch:

1. (Optional) Ensure your fork's `main` branch is up-to-date with the `AlexsLemonade/OpenScPCA-analysis` repository's main branch.
    - Refer to [this documentation about how to keep your `main` branch in sync](STUB_LINK stay in sync).

2. Open GitKraken and make sure you are on the `main` branch of the upstream repository, `AlexsLemonade/OpenScPCA-analysis`.
You can do this by double-clicking on the `main` branch in the branch graph with the ALSF logo.

<figure markdown="span">
    ![Click button to checkout main](../../img/working-with-branches-1.png){width="600"}
</figure>

3. Create and checkout a new branch by right-clicking on the `main` branch and selecting `Create branch here`.

<figure markdown="span">
    ![Create branch](../../img/working-with-branches-2.png){width="600"}
</figure>

4. You will then be prompted to name your branch.
We encourage users to name their branch using their GitHub username followed by a description of what changes will be included in the branch.
For example, if you are adding an analysis to add cell type annotations to Ewing sarcoma samples, you might name your branch `username/celltype-ewings`.

<figure markdown="span">
    ![Name branch](../../img/working-with-branches-3.png){width="600"}
</figure>

!!! note "Helpful hints for choosing a branch name"

    - Branch names should always be unique, so don't reuse old branch names!
    - Names should not contain any special symbols besides `_` and `-`.
    - It can be extra helpful to include the issue number in your branch name, e.g., `username/42-celltype-ewings`.

You have now successfully created your new feature branch.
You can now proceed to make changes to the code and continue with your analysis!
