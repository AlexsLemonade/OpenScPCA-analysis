
# Install WSL 2

!!! note
    These instructions apply to Windows users _only_.


## Prerequisites

WSL 2 requires that you are on Windows 11 or Windows 10 version 2004+ (Build 19041 and higher).
If you do not have at least these Windows versions, you must use [Lightsail for Research](../software-platforms/aws/index.md#lightsail-for-research-virtual-computing-with-aws) instead.


## Installation instructions

These instructions will provide you with both WSL 2 and the Ubuntu terminal:

- WSL 2 allows your computer to have a separate Linux subsystem
- Ubuntu is the default Linux distribution that gets installed when you install WSL 2
    - The Ubuntu terminal application also acts as the Linux [terminal interface](../getting-started/project-tools/using-the-terminal.md).
    - The Ubuntu terminal is the main way you will interact with the Linux side of your computer.
        - To learn more about the relationship between your Linux and Windows file systems, please refer to this article on [Working across Windows and Linux file systems](https://learn.microsoft.com/en-us/windows/wsl/filesystems).
        However, please be aware that for the purposes of OpenScPCA contribution, you should _only_ install files and perform analyses on the Linux file system.

To install WSL 2, take the following steps.
Note that these instructions are taken from [these official Windows instructions](https://learn.microsoft.com/en-us/windows/wsl/setup/environment).

Throughout this process you will be prompted about whether you want to all this app to make changes to your device.
Always click "Yes" when you see these prompts.

1. In the Windows Search Bar Menu, search for the "Windows PowerShell" application.
Open it by clicking "Run as administrator".

1. Run this command in PowerShell to install WSL 2:

    ```sh
    wsl --install
    ```

1. Once WSL 2 has finished installing, you will be prompted to set up a username and password to use with Ubuntu.
    - These credentials are _independent_ of the username and password you already have set up on the Windows side of your computer.
        - Changing one will not affect the other, but you can use the same username for both if you would like.
    - _Make sure you keep track of your username and password!_
    You will need to use your password when installing software for Ubuntu, and if you choose to use [RStudio](environment-setup/install-r-rstudio.md#using-the-rstudio-server), you will need your username and password.
    - Note that when you type your password, no symbols will appear - this is expected!

1. Open an [Ubuntu terminal window](../getting-started/project-tools/using-the-terminal.md#accessing-the-terminal-on-wsl2-on-windows).

    - Run the following command in Ubuntu to ensure that the package index for `apt`, the native Ubuntu package manager, and all its pre-installed packages are up to date.
    - Ubuntu will prompt you for your newly-created password when you run this command; you can expect to be prompted for a password any time you run a command as [`sudo`](https://www.pluralsight.com/resources/blog/cloud/linux-commands-for-beginners-sudo).
         - Again, when you type your password, no symbols will appear, as expected.
         - Enter "Y" if/when you are prompted for whether you want to continue with any package upgrades.

    ```sh
    sudo apt update && sudo apt upgrade
    ```


!!! tip "An Ubuntu tip!"
    Be aware that right-clicking in the Ubuntu terminal window has different behavior compared to the rest of your computer.
    In Ubuntu, right-clicking will _paste_ the contents of your clipboard into the terminal.
