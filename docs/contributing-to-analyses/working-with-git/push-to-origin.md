# Pushing to origin

Pushing commits from your local branch to a remote branch creates a remote copy of the changes that you made on your local branch.
This means your changes can now be viewed locally and on GitHub.

!!! note "Learn about pushing to origin"

    - [Pushing commits to a remote repository](https://docs.github.com/en/get-started/using-git/pushing-commits-to-a-remote-repository)
    - [Pushing to origin with GitKraken](https://www.gitkraken.com/learn/git/problems/git-push-to-remote-branch)

## Why push?

Creating a remote copy of your branch makes it possible to [file a pull request](STUB-LINK) to the `main` branch of `OpenScPCA-analysis`.
Filing a pull request is the only way to incorporate your changes into the main code base of OpenScPCA.


## How to push

1.Make sure all changes you would like to include have been committed.
See the section on [making commits](./making-commits.md).

2.Select `Push`, located above the branch diagram in `GitKraken`.

<figure markdown="span">
    ![Click push](../../img/push-to-origin-1.png){width="600"}
</figure>

The first time that you push to origin for a feature branch, you will be prompted to select the remote branch to push to and from.

- The first box should say `origin`, telling GitKraken to push your changes to a branch to your remote fork of `OpenScPCA-analysis`.
- The second box tells GitKraken what to name the remote copy of the branch.
    - This field will be automatically populated with the same name used for the local branch.
    - It is recommended to keep this and use the same branch name for the local and remote branch.
- Press `Submit` to confirm your choices and push to origin.

<figure markdown="span">
    ![Submit remote branch name](../../img/push-to-origin-2.png){width="600"}
</figure>

3.Confirm that you now have a remote copy of your feature branch.

After you push to origin, you will see a message pop-up in the lower left stating, `Pushed successfully`.

<figure markdown="span">
    ![Pushed successfully](../../img/push-to-origin-3.png){width="600"}
</figure>

You will also be able to see both the local computer icon and your GitHub avatar next to your branch name, indicating that your local branch and remote branch are now in sync.

<figure markdown="span">
    ![Sync local and remote branches](../../img/push-to-origin-4.png){width="600"}
</figure>

4.If more changes are needed, commit those changes, and then select `Push` to push to origin.

Because you have already set the remote branch, changes will automatically be pushed to that remote branch using the `Push` button, with no further prompts from GitKraken.

Congratulations, you have now synced your local changes with a remote branch!
