# Set up pre-commit

## What is pre-commit?

[`pre-commit`](https://pre-commit.com/) is a small software package that makes it easy to manage and run code quality checks prior to committing code.

As an OpenScPCA contributor, you will use pre-commit when you [commit files](../../contributing-to-analyses/working-with-git/making-commits.md) while writing analysis code.

`pre-commit` checks code quality by defining a set of "hooks," called pre-commit hooks, that will automatically run every time you commit changes to a repository.

- These pre-commit hooks check that the changes you are trying to commit pass the code quality standards we have set up.
- If `pre-commit` finds that certain files do not pass these checks, you will need to update the files to ensure they pass before you can commit your changes.

  - You can learn more about this process in the [contributing to analyses documentation](STUB_LINK).

## How does OpenScPCA use pre-commit?

All OpenScPCA contributors should use `pre-commit` while [making contributions with Git](./../../contributing-to-analyses/working-with-git/index.md).
Briefly, you will follow process:

- Make code changes [in a feature branch](../../contributing-to-analyses/working-with-git/working-with-branches.md)
- [Stage and Commit](../../contributing-to-analyses/working-with-git/making-commits.md) files with code changes
    - ‼️ Fix errors flagged by the pre-commit hook until the commit is successful
- [Push changes](../../contributing-to-analyses/working-with-git/push-to-origin.md) to your fork on GitHub
- [File a pull request](STUB_LINK for filing a PR) to the `AlexsLemonade/OpenScPCA-analysis` upstream repository


We have set up some pre-commit hooks to manage basic code security and catch other common problems, such as the following:

- Large data files that should not be committed to the repository
- [Merge conflicts](STUB_LINK resolving merge conflicts) that have not yet been resolved
- Credential files and other sensitive information


## Setting up pre-commit

During [conda setup](./setup-conda.md/#set-up-conda), you should have installed the `pre-commit` software into your base environment.

To turn on the pre-commit hooks for the OpenScPCA repository, you will need to run the command `pre-commit install` from a terminal window inside the repository.
We recommend launching [a terminal window in GitKraken](../../software-platforms/general-tools/using-the-terminal.md#how-do-you-access-the-terminal) to do this.

1. Open a terminal window and navigate to the `OpenScPCA` repository.
    - Remember, if you open the terminal within GitKraken, you're automatically in the repository - no extra navigation steps needed!


1. Type this command in the terminal prompt, and hit enter:

    ```sh
    pre-commit install
    ```

1. The terminal should then return this output message showing that you successfully set up your pre-commit hooks:

    ```
    pre-commit installed at .git/hooks/pre-commit
    ```

All set! 🎉


## Advanced considerations

_Instructions in this section are optional._

While we use a limited set of required pre-commit hooks in this project, there are [several other pre-commit hooks that you might find useful for your own development](#suggested-additional-pre-commit-hooks).
This section describes how you can add additional pre-commit hooks.

Because all OpenScPCA contributors share the `.pre-commit-config.yaml` file, please do not modify it.
Instead, you should make a copy of this file in the root directory of the repository named `.pre-commit-local.yaml`.
We have configured Git to ignore this file in the repository's [`.gitignore` file](https://docs.github.com/en/get-started/getting-started-with-git/ignoring-files), so you can modify it as you wish without affecting other contributors.

To switch your local activation of `pre-commit` to use this file, you will need to re-activate `pre-commit` using the the following code:

```sh
# make and activate a local pre-commit configuration
cp .pre-commit-config.yaml .pre-commit-local.yaml
pre-commit install --config .pre-commit-local.yaml
```

You can then add your own pre-commit hooks by editing the `.pre-commit-local.yaml` file.

### Suggested additional pre-commit hooks

#### Code formatting and linting

Code formatting and linting hooks can help you write error-free, consistent, and readable code.
For more on the value of these tools, [see this article about linters and formatters](https://www.freecodecamp.org/news/using-prettier-and-jslint/) (although this article focuses on JavaScript, the same principles apply to other languages).

The OpenScPCA required pre-commit hooks do not include any code formatting or linting checks, but you may find it helpful to use these tools in your own copy of the repository.

Note that these tools will often directly modify your files if they find errors.
In this case, the initial commit will fail, and you will then need to check and stage the changes that the tool made before retrying the commit.

Some formatters that we recommend are [`ruff-format`](https://docs.astral.sh/ruff/formatter/) for Python and [the `style-files` hook from the pre-commit package](https://lorenzwalthert.github.io/precommit/articles/available-hooks.html#style-files) for R.

You can add these hooks by adding the following code to your `.pre-commit-local.yaml` file in the `repos:` section:

```yaml
  # ruff formatter for Python
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.3.3
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

Code linting tools are often more intrusive compared to code formatting tools.
For example, they enforce not only general formatting but also particular style standards and "best" practices.
This can make them more likely to find errors and inconsistencies, and therefore also more likely to require you to manually fix the errors it finds.
They also might flag things that are not actually errors, but are simply not to the linter's taste.

If you would like to use a code linter, we recommend using [`ruff`](https://docs.astral.sh/ruff/) (which goes along with `ruff-format` shown above) for Python and the [`lintr` package](https://lintr.r-lib.org/) for R:

```yaml
  # ruff linter for Python
  - repo: https://github.com/astral-sh/ruff-pre-commit
    rev: v0.3.3
    hooks:
      - id: ruff
  # code linting with the  ith the {lintr} package for R
  - repo: https://github.com/lorenzwalthert/precommit
    rev: v0.4.0
    hooks:
      - id: lintr
```

#### Spell checking

Spell checking is another useful class of tools that can be added as a pre-commit hook.
Because biological and computational words are often not in default dictionaries, it is most helpful to use a spell checker that looks for common errors rather than one that checks every word against a dictionary.

One such tool is `typos`, which runs quickly to check both text and code for common mistakes.
`typos` can be installed as a pre-commit hook with the following code added to the `.pre-commit-local.yaml` file:

```yaml
  - repo: https://github.com/crate-ci/typos
    rev: v1.18.2
    hooks:
      - id: typos
```

#### Other pre-commit hooks

There are many other pre-commit hooks available.

You can find a fairly extensive list of them at the pre-commit ["Supported hooks" page](https://pre-commit.com/hooks.html).
Just note that the more hooks you add, the longer each commit to your repository will take!


#### Updating pre-commit hooks

If you ever need to update the versions of your pre-commit hooks, please use this command:

```sh
pre-commit autoupdate -c .pre-commit-local.yaml
```