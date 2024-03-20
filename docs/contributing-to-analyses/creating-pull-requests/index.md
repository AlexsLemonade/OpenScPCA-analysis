# Creating pull requests

## What is a pull request?

A pull request is used to propose merging new changes saved in one branch into another branch.
In most cases, pull requests are used to merge in changes from a [feature branch](../working-with-git/working-with-branches.md) into the `main` branch of the repository.

Pull requests can be viewed on GitHub, and include:

- A summary of line-by-line differences between the feature branch and the main branch.
- A summary, written by the contributor, describing the contents of changes proposed in the pull request.
- A history of all commits made on the feature branch.

!!! note "Learn more about pull requests"
    For more details on pull requests, see:

    - [GitHub's documentation describing pull requests](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/about-pull-requests)
    - [Tutorial on pull requests from GitKraken](https://www.youtube.com/watch?v=2VX1ISk9XH8&t=1s)

## Using pull requests in OpenScPCA

Once you have created a [feature branch](../working-with-git/working-with-branches.md) and [committed changes to that branch](../working-with-git/making-commits.md), you will file a pull request to merge those changes into the `main` branch of `AlexsLemonade/OpenScPCA-analysis`.

- Each pull request contains a group of related changes to the code.
    - Not sure what to include in your pull request?
    See [scoping a pull request](STUB-LINK to scoping).
- Each pull request includes a description of the changes made.
    - We provide a [pull request template](STUB-LINK to templates) to guide contributors on what to include in their pull request description.
- All pull requests are [reviewed](STUB-LINK to review) by at least one Data Lab staff member before they can be approved.
Once approved, changes are merged into the `main` branch of `AlexsLemonade/OpenScPCA-analysis`.

## Reviewing pull requests

We require that the content of all pull requests undergo [review by a Data Lab staff member](STUB-LINK to reviewing).
This includes:

- Checking for correctness, clarity, and reproducibility in the code
- Evaluating methods and rationale for any proposed analysis
- Ensuring that all necessary documentation is present and clear

If the reviewer has comments or requests changes, those comments should be addressed prior to requesting a second review.
Once the reviewer is happy with the changes, they will approve the pull request.

Approved pull requests can be safely [merged and incorporated](STUB-LINK merging).
Any new changes will be present in the `main` branch of `AlexsLemonade/OpenScPCA-analysis`.
