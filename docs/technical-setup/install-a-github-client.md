# Install a GitHub client

You will need a platform on your computer to run Git commands and interact with [GitHub.com](https://github.com).

There are many platforms you can choose, but we recommend using either [**GitKraken**](https://www.gitkraken.com/) or the [Git **command line interface**](https://git-scm.com/book/en/v2/Getting-Started-The-Command-Line).

## Why use GitKraken?

!!! note
    This option is best for contributors who are new to Git.
    The OpenScPCA documentation will primarily present how to use Git via GitKraken.

GitKraken is a GUI (graphical-user interface) which lets you "point-and-click" to run all Git commands.
GitKraken is a more user-friendly option which does not require you to use the command line (aka, terminal) or memorize Git commands.
It also provides a more streamlined approach to using Git when working on a Windows machine, v.s. using Gitcommand line interface.
GitKraken is also free to use with public repositories like `OpenScPCA-analysis`; you do not need to upgrade to a paid plan.

### Install GitKraken

1. Install GitKraken [using this link](https://www.gitkraken.com/download).

2. Set up GitKraken on your machine by [directly signing in with your GitHub account](https://help.gitkraken.com/gitkraken-client/github-gitkraken-client/#sign-in-with-github).
This will automatically provide you with the credentials you need to interact with GitHub without further setup.


## Why use the Git command line interface?

!!! note
    This option is best for contributors who have previous Git and command-line experience.
    The OpenScPCA documentation will generally not present how to use Git via the command line.<br><br>
    If you are new to Git or have limited command-line experience, we do not recommend this option.

The UNIX command-line interface is the classic way of running Git commands.
We specifically recommend this option for contributors who have previous experience working with Git on the command line.

If you choose this option, you will need to install Git itself if you do not already have it.
The instructions for this depend on your operating system, as described below.

After installing Git, you will also need take these steps:

- [Set up an SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh) to interact with GitHub
- [Set up your Git config file](https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup)

### Install Git on MacOS

We recommend installing Git via Apple's developer command line tools.
For this, launch a Terminal window, enter the command shown below, and follow the instructions.

```
xcode-select --install
```


### Install Git on Windows

We recommend installing the [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install) for a full UNIX command-line experience via Ubuntu.
You can then [install Git](https://git-scm.com/download/linux) into your Ubuntu system.


