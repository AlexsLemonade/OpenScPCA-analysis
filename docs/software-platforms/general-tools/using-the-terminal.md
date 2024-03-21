# Using the terminal

This page provides some background information on working with the terminal, also known as the command line prompt.


## What is the terminal?

The terminal is an interface where you can type commands (rather than "point and click") to interact with your operating system.
In the OpenScPCA project, we use the word "terminal" to specifically refer to a UNIX-style terminal, and not for example the Windows PowerShell.

macOS and Linux (which [Lightsail for Research instances](STUB_LINK to LSfR docs) use) computers come with a built-in terminal application.
On Windows computers, you will first need to install and setup the [Windows Subsystem for Linux (WSL)](STUB_LINK for WSL instructions) to be able to use the terminal.

!!! info "Learn more about the terminal"

    If you are new to the terminal, we highly recommend starting with a general introduction.
    You do not need to become a terminal expert, but it's helpful to have some basic bearings.

    Here are a few we resources we suggest.
    Note that while some links have details that are specific to macOS or Linux, the general principles and commands are exactly the same on both kinds of computers!

    - [Introduction to the command line](https://tutorials.codebar.io/command-line/introduction/tutorial.html)
    - macOS
        - [How to Use Terminal on a Mac: A Beginner's Guide](https://www.makeuseof.com/tag/beginners-guide-mac-terminal/)
        - [Video: Using the Terminal on macOS for beginners](https://www.youtube.com/watch?v=dWFMRp6KtlQ)
    - Linux and WSL on Windows
        - [Command line for beginners](https://ubuntu.com/tutorials/command-line-for-beginners#1-overview)

## When will you need to use the terminal?

There are certain situations that you will need to use a terminal as part of contributing to OpenScPCA, such as:

- To perform the [initial conda setup](../../technical-setup/environment-setup/setup-conda.md)
- To set up the [`pre-commit` package](../../technical-setup/environment-setup/setup-precommit.md)
- To [create a new analysis module](../../contributing-to-analyses/analysis-modules/STUB_LINK create a module)
- To [run all full analysis module](STUB_LINK run an analysis module)

## How do you access the terminal?

First, no matter which kind of computer (Mac, WSL in Windows, or Lightsail for Research instance) you are using, you can directly access a terminal within GitKraken.
The benefit of using the GitKraken terminal is that it always places you in your `OpenScPCA-analysis` repository folder, saving you time from having to navigate to the repository in terminal with the `cd` command.

To launch a terminal in GitKraken, first open your repository in GitKraken.

Then, click the `Terminal` button at the top of GitKraken.
This will open a terminal prompt inside GitKraken where you can type your commands.

<figure markdown="span">
    ![Open the terminal in GitKraken.](../../img/gitkraken-terminal-button.png){width="600"}
</figure>


Or, you can launch the terminal application directly, as explained below for different operating systems.


### Accessing the terminal on a Mac

The macOS Terminal application is stored in `Applications -> Utilities`.

<figure markdown="span">
    ![Open the terminal in macOS.](../../img/launch-terminal-mac.png){width="600"}
</figure>

Double-click this application to open the terminal.
(PS: You can also change the look of the terminal by [customizing your terminal profile](https://support.apple.com/guide/terminal/profiles-change-terminal-windows-trml107/mac).)


### Accessing the terminal on WSL on Windows

Make sure that you have taken all the necessary steps to [install WSL](STUB_LINK WSL instructions):

- First, install WSL itself
- Second, install the Ubuntu app from the Microsoft Store

This Ubuntu app itself _is_ the terminal!
To launch the terminal, simply [search for and open the Ubuntu app via the Windows Search Bar Menu](https://canonical-ubuntu-wsl.readthedocs-hosted.com/en/latest/guides/install-ubuntu-wsl2/#method-1-microsoft-store).

### Accessing the terminal on Lightsail for Research

Lightsail for Research instances run on the Ubuntu operating system, which is a type of Linux.

To launch the terminal in Ubuntu, INSTRUCTIONS ARE FORTHCOMING FOR WHERE EXACTLY TO CLICK.
