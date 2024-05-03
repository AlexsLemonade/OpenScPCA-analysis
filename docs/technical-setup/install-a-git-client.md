# Install a Git client

You will need a platform on your computer to run Git commands and interact with [GitHub](https://github.com).

There are many platforms you can choose, but we recommend using either [**GitKraken**](https://www.gitkraken.com/) or the [Git **command line interface**](https://git-scm.com/book/en/v2/Getting-Started-The-Command-Line).

## Why use GitKraken?

!!! note
    This option is best for contributors who are new to Git.
    The OpenScPCA documentation will primarily present how to use Git via GitKraken.

GitKraken is a GUI (graphical-user interface) for using Git, meaning you can "point-and-click" to run all Git commands.
It is therefore a user-friendly option which does not require you to use the command line (aka, terminal) or memorize Git commands.

GitKraken is _free_ to use with public repositories like `OpenScPCA-analysis`; you do not need to upgrade to their paid plan!

Be aware that installing GitKraken will install Git, but _only_ for use in GitKraken.
If you want to also be able to use Git via command line, please refer to the [command line installation instructions](#why-use-the-git-command-line-interface).

### Install GitKraken on macOS

1. Install GitKraken [using this link](https://www.gitkraken.com/download).
Note that installing GitKraken will also provide you with Git itself.

1. Set up GitKraken on your machine by [directly signing in with your GitHub account](https://help.gitkraken.com/gitkraken-client/github-gitkraken-client/#sign-in-with-github).
This will automatically provide you with the credentials you need to interact with GitHub without further setup.

### Install GitKraken on Windows with WSL 2

??? note "Update WSL 2 if needed"
    Before installing GitKraken, make sure your WSL 2 installation is up to date by running this command in Windows PowerShell (not in Ubuntu).
    When launching PowerShell, make sure you are [running as an administrator](https://learn.microsoft.com/en-us/powershell/scripting/windows-powershell/starting-windows-powershell?view=powershell-7.4#with-administrative-privileges-run-as-administrator).

    ```sh
    wsl --update
    ```

    Still need to install WSL 2?
    [Follow these instructions first](./install-wsl-2.md).

You will [install GitKraken as a Linux application](https://help.gitkraken.com/gitkraken-client/how-to-install/#deb) into the WSL 2 side of your computer.

All of the following installation steps are commands that you should run in the Ubuntu terminal.


1. Download the GitKraken Linux installer:

    ```sh
    wget https://release.gitkraken.com/linux/gitkraken-amd64.deb
    ```

1. Install GitKraken:

    ```sh
    sudo apt install ./gitkraken-amd64.deb
    ```

    Along the way, you may be prompted `Do you want to continue [Y/n]?`.
    When you see this prompt, hit enter (or "Y") to continue with the installation.

1. You can now safely remove the installer file:

    ```sh
    rm gitkraken-amd64.deb
    ```

1. Now that you have installed GitKraken, you can open it by typing `gitkraken` in the Ubuntu terminal.

1. The final step is to set up GitKraken on your machine by [directly signing in with your GitHub account](https://help.gitkraken.com/gitkraken-client/github-gitkraken-client/#sign-in-with-github).
This will automatically provide you with the credentials you need to interact with GitHub without further setup.
    - Note that you may need to manually copy/paste the Token displayed in your browser into GitKraken to complete the sign-in.
    - You can also refer to the [GitKraken website](https://help.gitkraken.com/gitkraken-client/windows-subsystem-for-linux) for additional considerations (some of which are more advanced) for using GitKraken with WSL 2.



## Why use the Git command line interface?

!!! note
    This option is best for contributors who have previous experience with Git and the command line.
    The OpenScPCA documentation will generally not present how to use Git via the command line.<br><br>
    If you are new to Git, we do not recommend this option.

The command line interface is the classic way of interacting with Git.
We specifically recommend this option for contributors who have previous experience working with Git on the command line.

If you choose this option, you will need to install Git itself if you do not already have it.
The instructions for this depend on your operating system, as described below.

After installing Git, you will also need take these steps:

- [Set up an SSH key](https://docs.github.com/en/authentication/connecting-to-github-with-ssh) to interact with GitHub
- [Set up your Git config file](https://git-scm.com/book/en/v2/Getting-Started-First-Time-Git-Setup)

### Install Git on macOS

We recommend installing Git via Apple's developer command line tools.
For this, launch a Terminal window, enter the command shown below, and follow the instructions.

```
xcode-select --install
```

### Install Git on Windows in WSL 2

Git is automatically installed as part of WSL 2, so you don't need to install anything else.
