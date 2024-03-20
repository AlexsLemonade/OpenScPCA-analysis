# Scoping a pull request

## What to include in a pull request

When filing a PR, the content should be focused and contain a set of related changes.
For most cases, each PR should contain changes that address one [issue filed in the repository](../../communications-tools/github-issues/what-are-github-issues-and-how-do-we-use-them.md).

Keeping PRs focused on one task can help by combatting reviewer fatigue and keeping the [review process](STUB-LINK to review) short.
Longer PRs are harder to review and it can be much easier to introduce an unwanted bug in the code.
The smaller the PR, the faster and more thorough the review!

For example, if you are working on adding an [analyses module](STUB-LINK to analysis modules), you may file the following PRs:

- The first PR should include [initiating the analysis module](STUB-LINK to creating analysis module), including the skeleton of the module directory and a description in the `README.md` file.
- The second PR may include a script (and any accompanying documentation) that performs the first step of the analysis.
- Subsequent PRs can add in any additional scripts, where each PR adds a single step in the analysis.

## Rules of thumb for pull requests

Below are some good rules of thumb to follow when creating a PR:

- The total number of lines that are changed should be less than 400.

    When preparing to file a PR on GitHub, you can scroll below the template description to see the summary of commits and changes.
    Below the summary of all commits, a summary of your line-by-line changes with the number of lines with additions and deletions will be displayed.
    You can view this before filing your PR to ensure you have a PR that is appropriate length.

    <figure markdown="span">
        ![Commit summary](../../img/scoping-prs-1.png){width="600"}
    </figure>

- The contents of the PR address the contents of a single issue and solve one problem.
This means for each issue, there is one PR that addresses that issue.
_Please do not_ file PRs that address multiple issues at once.

