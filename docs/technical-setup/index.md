# Technical setup

This section provides instructions for technical aspects of set up, including:

- Downloading and setting up a [Git client](./install-a-git-client.md)
- [Forking the `AlexsLemonade/OpenScPCA-analysis` repository](./fork-the-repo.md)
- [Cloning your fork](./clone-the-repo.md) to your computer
- Setting up [additional dependencies](environment-setup/index.md) on your computer that you'll need to contribute to OpenScPCA

## Special considerations for Windows users

Due to the OpenScPCA project's technical environment needs, Windows users will [need to install and use Windows Subsystem for Linux 2 (WSL 2)](./install-wsl-2.md).


[WSL 2](https://learn.microsoft.com/en-us/windows/wsl/about) lets you run a separate Linux (Ubuntu) environment on your Windows computer.
This includes a separate Linux file system that is mostly distinct from your regular Windows file system.
When contributing to OpenScPCA, all data analysis will happen on the WSL 2 side of your computer, which means:

- Your repository needs to be cloned into your WSL 2 file system
- All applications (e.g., GitKraken and RStudio) need to be installed as Linux applications on your WSL 2 system

_The OpenScPCA project does not support the use of Windows machines without WSL 2._
Note that all contributors have access to [Lightsail for Research (LSfR)](../software-platforms/aws/index.md#lsfr-virtual-computing-with-aws) instances which run Linux (specifically, Ubuntu), so you can use LSfR if you prefer not to install WSL 2.

<!-- TODO: We have also created a specific GitHub Discussions category where you can post questions about using WSL 2.-->

!!! note "Learn more about WSL 2"
    For more information about WSL 2, please see this [Introduction to Windows Subsystem for Linux](https://learn.microsoft.com/en-us/training/modules/wsl-introduction/).

