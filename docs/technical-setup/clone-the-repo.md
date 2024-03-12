# Clone the repository

After [forking the `OpenScPCA` repository](./fork-the-repo.md), you will need to _clone_ your forked repository to your computer.
Cloning provides you with a local copy of the repository to work with.
You will need to clone the repository to every computer (or remote server) you want to use to contribute to the project.

!!! note
    You should clone your forked repository from `<username>/OpenScPCA-analysis`, _not_ `AlexsLemonade/OpenScPCA-analysis`.


## Clone your forked repository

Please follow the [instructions in this GitKraken video](https://help.gitkraken.com/gitkraken-client/open-clone-init/#cloning-an-existing-project) to clone your repository.

This video presents two approaches to cloning: `Clone with URL`, and Clone based on your integrated GitHub account.
Either approach is fine to take, but please bear in mind the following:

In the `Clone with URL` option, the video presents cloning with an `https` link, vs. an `SSH` link.
If you have already set up a [GitHub SSH connection](https://docs.github.com/en/authentication/connecting-to-github-with-ssh), you are welcome to use the `SSH` link version instead.
If you do not have SSH set up for GitHub, you will want to stay with the `https` link as shown.
Note that an SSH connection is _not required_ to contribute to OpenScPCA.

While watching this video, you can feel free to ignore all references to `LFS` ("large file storage"); this does not apply to `OpenScPCA-analysis`.

## Add the project repository as a remote repository

The next step is to link the upstream repository (`AlexsLemonade/OpenScPCA-analysis`) as a remote.
A _remote_ is a repository on GitHub that you are connected to.
Each remote repository is given a name to make referring to them easier.
For example, your forked repository on GitHub is called `origin` because it is where your local repository was cloned from.

We will add another remote repository named `upstream` that refers to the original `AlexsLemonade/OpenScPCA-analysis` repository on GitHub.
This will be called `upstream` because it was the source that your fork came from.

Adding the `AlexsLemonade` _upstream remote_ will allow you to interact with it from your computer which can help you keep your fork in sync with the `OpenScPCA` project.
But, you will still be working in your fork when writing analysis code.

Follow these steps to add the upstream remote:

1. From your repository in GitKraken, hover over the `1/1` text on the left-hand side `Remote` menu.
This text will then turn into a plus-sign icon.
Click that icon.
    <!-- keep this tabbed in to enable the numbered list -->
    <figure markdown="span">
        ![Click button to add the remote.](../img/add-upstream-remote-1.png){width="600"}
    </figure>

1. The following screen will prompt you to add a remote.
Select `AlexsLemonade/OpenScPCA-analysis` from the dropdown menu, and click the button `Add remote`.
    <figure markdown="span">
        ![Add the upstream remote.](../img/add-upstream-remote-2.png){width="400"}
    </figure>

1. You should then see a second remote called `AlexsLemonade` on the left-hand side `Remote` panel.
The specific listed items under `AlexsLemonade` that you see will look different from the screenshot below; this is expected.
As long as you see that `AlexsLemonade` is listed in the menu, you have successfully added the upstream remote.
    <figure markdown="span">
        ![View the added remote.](../img/add-upstream-remote-3.png){width="400"}
    </figure>
