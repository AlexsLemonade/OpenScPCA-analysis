# Staying in sync with upstream branches

You may wish to keep changes made in `AlexsLemonade/OpenScPCA-analysis` in sync with the `main` branch of your fork.
This ensures that your fork has the most recent changes and updates found in the main code base.

!!! note "Learn about syncing your fork"

    For more information on syncing your fork with the upstream repository, see:

    - [GitHub documentation on syncing a fork](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/syncing-a-fork)


## How to sync your fork with changes in `AlexsLemonade/OpenScPCA-analysis`

Follow these steps to keep your fork in sync with changes in `AlexsLemonade/OpenScPCA-analysis`:

1. Login to GitHub and navigate to the main page of your forked repository.

1. If the `main` branch of your fork is out of date with the `main` branch of `AlexsLemonade/OpenScPCA-analysis`, then you will see a note that states, "This branch is X commits behind `AlexsLemonade/Open-ScPca-analysis:main`".

    You can then press `Sync fork` to sync your fork with any upstream changes.

    <figure markdown="span">
        ![Sync fork](../../img/upstream-sync-1.png){width="600"}
    </figure>

1. After pressing `Sync fork` a pop-up will notify you that the branch is out-of-date.
    Go ahead and press `Update branch`.

    <figure markdown="span">
        ![Update branch](../../img/upstream-sync-2.png){width="400"}
    </figure>

1. If your `main` branch has been updated successfully, you will now see a note that the branch is up to date with `AlexsLemonade/OpenScPCA-analysis`.

    <figure markdown="span">
        ![Successful sync](../../img/upstream-sync-3.png){width="600"}
    </figure>

## Resolving sync conflicts

If there are changes that you have on your `main` branch that conflict with changes made in the `main` branch of `AlexsLemonade/OpenScPCA-analysis`, you will not be able to sync until those conflicts have been resolved.

When this happens, GitHub will notify you that there are conflicts after pushing `Sync fork`.

<figure markdown="span">
    ![Sync conflicts](../../img/upstream-sync-4.png){width="600"}
</figure>

Follow these steps to manage sync conflicts with GitKraken:

1. Open up your forked repository in GitKraken and checkout the `main` branch of your fork.

1. Repeat the attempt to sync your changes with the upstream branch.
    To do this, right-click on the remote `main` from `AlexsLemonade/OpenScPCA-analysis` and select `merge AlexsLemonade/main into main`.

    <figure markdown="span">
        ![Merge main](../../img/upstream-sync-5.png){width="400"}
    </figure>

1. A banner will pop up indicating that there is a merge conflict and any files that contain conflicts will be listed on the right-hand side under `Conflicted files`.

    <figure markdown="span">
        ![Conflict banner](../../img/upstream-sync-6.png){width="600"}
    </figure>

    <figure markdown="span">
        ![Conflict files list](../../img/upstream-sync-7.png){width="600"}
    </figure>

    For help in resolving these conflicts, please watch this [tutorial from GitKraken](https://www.gitkraken.com/learn/git/tutorials/how-to-resolve-merge-conflict-in-git).

    Once you have resolved the conflicts and committed the changes, your `main` branch will now be in sync with the upstream `main` branch of `AlexsLemonade/OpenScPCA-analysis`.
