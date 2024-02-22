# Contributing to OpenScPCA

**Table of Contents**
- [Setting up pre-commit](#setting-up-pre-commit)
  - [Additional optional hooks](#additional-optional-hooks)
    - [Code formatting and linting](#code-formatting-and-linting)
    - [Spell checking](#spell-checking)
    - [Other pre-commit hooks](#other-pre-commit-hooks)

## Setting up pre-commit

[`pre-commit`](https://pre-commit.com) is a small software package that makes it easy to manage and run code quality checks.
All contributors should use pre-commit as part of their workflow, installing the package as described below.
`pre-commit` checks code quality by defining a set of "hooks" that will run every time you commit changes to a repository.
We have used it in this project to set up some pre-commit hooks to manage basic code security and other common errors, such as the following:
- Large data files that should not be committed to the repository
- Credential files and other sensitive information
- Merge conflicts that have not yet been resolved

To install pre-commit, follow the instructions in the [pre-commit documentation](https://pre-commit.com/#install).
Once the `pre-commit` software is installed, copy our base pre-commit configuration file to `.pre-commit-config.yaml` and
activate the hooks by running the following commands in the root directory of the repository:

```bash
# run in the root directory of the repository
cp pre-commit-base.yaml .pre-commit-config.yaml
pre-commit install
```

After that point, any time you commit a change to the repository, the hooks will run and check for errors.
If any errors are found, the commit will be aborted and you will be prompted to fix the errors, after which you can retry the commit.
For some hooks, the errors will be automatically fixed, and you will only need to stage the updated files and retry the commit.

Note that the first time you commit after installing `pre-commit` or updating the hooks, the hooks may take a while to complete, as software may need to be downloaded and installed, but the hooks will run much faster on subsequent commits.

### Additional optional hooks

While we have taken a limited approach to the required pre-commit hooks in this project, there are a number of other pre-commit hooks that you might find useful for your own development.

You can add your own pre-commit hooks by editing the `.pre-commit-config.yaml` file in the root directory of the repository.
Changes to the `.pre-commit-config.yaml` file will be ignored by Git based on our `.gitignore` file, so you are free to modify it as you wish without affecting other users.

#### Code formatting and linting

One example is code formatting and linting tools, which can help you write error-free, consistent, and readable code.
For more on the value of these tools, see [this article about linters and formatters](https://www.freecodecamp.org/news/using-prettier-and-jslint/).
While the article focuses on JavaScript, the same principles apply to other languages.

We have not included any code formatting or linting checks in the pre-commit hooks we require, but you may find it helpful to use these tools in your own copy of the repository.

Note that these tools will often directly modify your files when run.
If they are run as a pre-commit hook the initial commit will fail, and you will then need to check and stage the changes that were made by the tool before re-trying the commit.


Some formatters that we recommend are [`ruff-format`](https://docs.astral.sh/ruff/formatter/) for Python and the [`style-files` hook from the precommit package](https://lorenzwalthert.github.io/precommit/articles/available-hooks.html#style-files) for R.
You can add those with the following code added to the `.pre-commit-config.yaml` file:

```yaml
  # ruff formatter for Python
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.2.1
    hooks:
      - id: ruff-format
  # code styling with the {styler} package for R
  - repo: https://github.com/lorenzwalthert/precommit
    rev: v0.4.0
    hooks:
      - id: style-files
  # prettier formatter for many other languages
  - repo: https://github.com/pre-commit/mirrors-prettier
    rev: v3.1.0
    hooks:
      - id: prettier
```

Code linting tools are often more intrusive, enforcing not only general formatting but also particular style standards and "best" practices.
This can make them more likely to find errors and inconsistencies, but also more likely to require manual intervention to fix those errors.
They also might complain about things that are not actually errors, but are simply not to the linter's taste.
For python, we recommend [`ruff`](https://docs.astral.sh/ruff/) (which goes along with `ruff-format`, above), and for R we recommend the [`lintr`](https://lintr.r-lib.org) package.


#### Spell checking

Spell checking is another useful class of tools that can be added as a pre-commit hook.
Because biological and computational words are often not in default dictionaries, it may be more helpful in a pre-commit context to use a spell checker that looks for common errors rather than one that checks every word against a dictionary.
One such tool is [`typos`](https://github.com/crate-ci/typos), which runs quickly, checking both text and code for common mistakes.
 `typos` can be installed as a pre-commit hook with the following code added to the `.pre-commit-config.yaml` file:

```yaml
  - repo: https://github.com/crate-ci/typos
    rev: v1.18.2
    hooks:
      - id: typos
```

#### Other pre-commit hooks

There are many other pre-commit hooks available, and you can find a fairly extensive list of them at the pre-commit ["Supported hooks" page](https://pre-commit.com/hooks.html).
Just note that the more hooks you add, the longer each commit to your repository will take!
