# Troubleshooting commits

In some cases, GitKraken will not allow you to make a commit.
Typically this will be because pre-commit hooks have not passed.

!!!note "More information on pre-commit"

    We ask that all contributors set up pre-commit hooks prior to making any commits.
    For more information see the [documentation provided on setting up pre-commit](../../technical-setup/setup-precommit.md).

## Why is pre-commit failing?

We have set up pre-commit hooks to manage basic code security and catch other common problems, such as:

- Large data files that should not be committed to the repository (files > 200 Kb)
- Merge conflicts that have not yet been resolved

    - If you need help resolving merge conflicts, please see [resolving merge conflicts](STUB-LINK for merge conflicts).

- Credential files and other sensitive information

If you attempt to commit anything in the above list, then you will see a red banner indicating pre-commit failed.
If this is the case, you will not be able to commit your changes until you fix the conflict.

You can click the `View log` button on the right hand side of the pre-commit banner to see which checks failed.

<figure markdown="span">
    ![Pre-commit fail](../../img/troubleshooting-commits-1.png){width="600"}
</figure>

In the below example, the pre-commit failure is because we are trying to commit a file that contains credentials.

<figure markdown="span">
    ![Pre-commit log](../../img/troubleshooting-commits-2.png){width="600"}
</figure>

The log may also include some details on which file caused an error and the reason.
Typically, this will require you to delete sensitive information or remove any large files.

Once you have identified and resolved the issue, you can commit your files and pre-commit should now pass.

<figure markdown="span">
    ![Pre-commit pass](../../img/making-commits-5.png){width="600"}
</figure>


## I added a file but it isn't showing up in GitKraken for me to stage and commit?

In some cases, you may go to commit a file that is not showing up in the `Unstaged files` section of GitKraken.
If the file is not there, than you cannot stage it and make a commit.

We have set up a `.gitignore` file in the root directory of this repository that will tell Git to ignore certain files.
Any files with a matching extension or path listed in the `.gitignore` will not be allowed to be commited to the repository.

This includes:

- Data files, such as `.rds` or `.hdf5` files
- Environment files, such as `renv` and `env` directories
- Other hidden files produced by R and Python
