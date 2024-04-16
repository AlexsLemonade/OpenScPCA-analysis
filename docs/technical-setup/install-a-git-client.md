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

### Install GitKraken on macOS

1. Install GitKraken [using this link](https://www.gitkraken.com/download).
Note that installing GitKraken will also provide you with Git itself.

1. Set up GitKraken on your machine by [directly signing in with your GitHub account](https://help.gitkraken.com/gitkraken-client/github-gitkraken-client/#sign-in-with-github).
This will automatically provide you with the credentials you need to interact with GitHub without further setup.

### Install GitKraken on Windows with WSL2

!!! note
    Before installing GitKraken, make sure your WSL2 installation is up to date by running this command in Windows Powershell (not in Ubuntu).
    When launching Powershell, make sure you are [running as an administrator](https://learn.microsoft.com/en-us/powershell/scripting/windows-powershell/starting-windows-powershell?view=powershell-7.4#with-administrative-privileges-run-as-administrator).

    Still need to install WSL2?
    [Follow these instructions first](./install-wsl2.md).

You will [install GitKraken as a Linux application](https://help.gitkraken.com/gitkraken-client/how-to-install/#deb) into the WSL2 side of your computer.
Note that installing GitKraken will also provide you with Git itself.

All of the following installation steps are commands that you should run in the Ubuntu app.

1. To begin, open the Ubuntu app and run this command:

    ```sh
    sudo apt update
    ```

    This command ensures that `apt`, the native Ubuntu package manager, is up to date.
    Note that Ubuntu will prompt you for your WSL2 password; again, no symbols will appear when you type it.

1. Download the GitKraken Linux installer:

    ```sh
    wget https://release.gitkraken.com/linux/gitkraken-amd64.deb
    ```

1. Install GitKraken:

    ```sh
    wget sudo apt install ./gitkraken-amd64.deb
    ```

    Along the way, you may be prompted `Do you want to continue [Y/n]?`.
    When you see this prompt, hit enter (or "Y") to continue with the installation.

Now that you have installed GitKraken, you can open it by typing `gitkraken` in the Ubuntu terminal.

The final step is to set up GitKraken on your machine by [directly signing in with your GitHub account](https://help.gitkraken.com/gitkraken-client/github-gitkraken-client/#sign-in-with-github).
This will automatically provide you with the credentials you need to interact with GitHub without further setup.

You can also refer to the [GitKraken website](https://help.gitkraken.com/gitkraken-client/windows-subsystem-for-linux) for additional considerations for using GitKraken with WSL2.



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

### Install Git on Windows in WSL2

Once you have [installed WSL2](./install-wsl2.md), you can install Git using `apt`:

```sh
# First make sure apt is up to date
sudo apt update

# Next, install Git
apt install git
```
