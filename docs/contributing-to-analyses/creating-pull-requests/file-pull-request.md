# Filing a pull request

Once you have [committed changes](../working-with-git/making-commits.md) to your [feature branch](../working-with-git/working-with-branches.md) and [pushed to origin](../working-with-git/push-to-origin.md), you are ready to [file a PR](index.md).

Don't know what to include in your PR or when to file a PR?

- See our [guidelines on scoping a PR](./scoping-pull-requests.md).

!!! note "More information on filing a pull request"

    For more on filing a pull request, see the [GitHub documentation on filing a PR](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request).

## How to file a pull request for OpenScPCA

All PRs should be filed from your feature branch to the `main` branch of `AlexsLemonade/OpenScPCA-analysis`.
Filing a PR is the only way to incorporate any new analysis or code changes into the main code base of OpenScPCA.

Follow the below steps to file a pull request on GitHub:

1. [Commit](../working-with-git/making-commits.md) all changes you plan to include in your PR to your [feature branch](../working-with-git/working-with-branches.md).

1. [Push](../working-with-git/push-to-origin.md) all commits to origin, ensuring that all commits are saved to the remote copy of your feature branch.

1. Navigate to your fork of the OpenScPCA repository on `GitHub.com`.

    Once there, use the drop-down menu on the top left to select the feature branch you would like to file the PR _from_.
    This will re-load the page and now display the contents of your feature branch.

    <figure markdown="span">
        ![Change branch](../../img/file-pr-1.png){width="300"}
    </figure>

1. Click on the `Contribute` button and select `Open pull request`.
Any changes found in your feature branch that are _not_ in the `main` branch of `AlexsLemonade/OpenScPCA-analysis` will be included in the PR.

    <figure markdown="span">
        ![Contribute](../../img/file-pr-2.png){width="300"}
    </figure>

1. You will then be redirected to a page where you can complete filing your PR.

    - This should display the base branch as the `main` branch of `AlexsLemonade/OpenScPCA-analysis`.
        - This is the branch you are proposing to merge your changes _into_.
    - The head repository should be your fork (`username/OpenScPCA-analysis`) with the branch set to your feature branch.
        - This is the branch containing the changes you would like to merge into the main code base of OpenScPCA.

    <figure markdown="span">
        ![Check branches](../../img/file-pr-3.png){width="600"}
    </figure>

    If you do not see a green checkmark with the words `Able to merge`, but instead see a red X and a warning stating `Merge conflicts detected`, you will need to [resolve merge conflicts](STUB-LINK for merge conflicts) before you can file the PR.

1. Be sure to provide a descriptive title summarizing your changes and fill out the [PR template](./pull-request-template.md).

    Please follow the instructions in the template and fill it out to completion.

    <figure markdown="span">
        ![Title and template](../../img/file-pr-4.png){width="600"}
    </figure>

1. Once you have completed the PR template, you can press `Create pull request`.

    After filing your PR, the Data Lab team will [review the proposed changes](STUB-LINK for review) before changes can be merged into the main code base.

You have now successfully created a pull request!
