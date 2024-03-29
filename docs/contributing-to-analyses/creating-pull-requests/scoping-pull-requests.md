# Scoping a pull request

When working on an analysis, it's important to take some time to plan your work.
For every new addition, you will need to [file a PR](./file-pull-request.md) and go through the [review process](../pr-review-and-merge/index.md).

- To ensure a fast and smooth review process, it is helpful to break down your proposed changes into small chunks or units of work.
    - We encourage users to [file issues](../../communications-tools/github-issues/index.md) describing each of these units of work!
- Each PR should contain changes related to a single unit of work.

For example, if you are working on adding an [analysis module](../analysis-modules/index.md), you may file the following PRs:

- The first PR should include [your initialized analysis module](../analysis-modules/creating-a-module.md), including the skeleton of the module folder and a description in the `README.md` file.
- The second PR may include a script (and any accompanying documentation) that performs the first step of the analysis.
- Subsequent PRs can add in any additional scripts, where each PR adds a single step in the analysis.

## Why are focused pull requests better?

All changes in a given PR should be focused with a set of related changes.
For most cases, each PR should contain changes that address one [issue filed in the repository](../../communications-tools/github-issues/index.md).
There may be occasions where one issue requires multiple PRs, but _do not_ file PRs that address multiple issues at once.

Ensuring your PRs are focused on one task helps combat reviewer fatigue and keeps the [review process](../pr-review-and-merge/index.md) short. <!--STUB_LINK: Replace with review process link -->

- Longer PRs are harder and more time-consuming to review, and reviewers are more likely to miss catching unwanted bugs in the code.
- The smaller the PR, the faster and more thorough the review!
- Quicker review means you are less likely to introduce [merge conflicts](resolve-merge-conflicts.md).

## Rules of thumb for good pull requests

Below are some good rules of thumb to follow when determining when to file a PR:

- The total number of changed lines should be less than 400.

    Any time that you are making changes on a branch, you can go through the steps of [filing a PR](./file-pull-request.md) to view the amount of lines that have been changed.
    Instead of actually filing the PR, you can scroll below the template description to see the summary of commits and changes.
    Below the summary of all commits, GitHub displays a summary of your line-by-line changes, including the number of lines with additions and deletions.

    In the example below, 18 lines have been changed.
    If you start seeing the number of lines approach 400, it's time to file that PR!

    <figure markdown="span">
        ![Commit summary](../../img/scoping-prs-1.png){width="600"}
    </figure>

- Each PR should address a single issue and solve one problem.
This means you should file, at a minimum, one PR to address each issue.
_Please do not_ file PRs that address multiple issues at once.

