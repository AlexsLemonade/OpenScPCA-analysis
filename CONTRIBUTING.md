# Contributing to OpenScPCA

**Table of Contents**
- [Pre-commit hooks](#pre-commit-hooks)
  - [Additional optional hooks](#additional-optional-hooks)

## Pre-commit hooks

We have set up pre-commit hooks to manage basic code security and other common errors using the [pre-commit framework](https://pre-commit.com).
We **strongly** encourage all contributors to use the included hooks as part of their workflow.
The included hooks check for things like the following:

To install pre-commit, follow the instructions in the [pre-commit documentation](https://pre-commit.com/#install).
Once pre-commit is installed, you can install the hooks by running `pre-commit install` in the root directory of the repository.

After that point, any time you commit a change to the repository, the hooks will run and check for errors.
If any errors are found, the commit will be aborted and you will be prompted to fix the errors, after which you can retry the commit.

Note that thie first time you commit after installing `pre-commit` or updating the hooks, the hooks may take a while to complete, as software may need to be downloaded and installed, but the hooks will run much faster on subsequent commits.

### Additional optional hooks

We have not included any code formatting or linting checks in the pre-commit hooks, but we encourage you to use these tools as well to ensure that your code is consistent and readable.
You can add to the pre-commit hooks by editing the `.pre-commit-config.yaml` file to the root directory of the repository.
Some formatters that we recommend are [`ruff-format`](https://docs.astral.sh/ruff/formatter/) for Python and the [`style-files` hook from the precommit package](https://lorenzwalthert.github.io/precommit/articles/available-hooks.html#style-files) for R.

You can add those with the following code added to the `.pre-commit-config.yaml` file:

```yaml
# ruff formatter for Python
- repo: https://github.com/astral-sh/ruff-pre-commit
  rev: v0.1.13
  hooks:
    - id: ruff-format
# code styling with the {styler} package for R
- repo: https://github.com/lorenzwalthert/precommit
  rev: v0.3.2
  hooks:
    - id: style-files
```

Note that these packages will modify the files in place when run, so you may need to add the changes they made to your commit after running the hooks.

Changes to the `.pre-commit-config.yaml` file will be ignored by Git based on our `.gitignore`` file, so you are free to modify it as you wish without affecting other users.

