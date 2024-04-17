
# Install WSL2

!!! note
    These instructions apply to Windows users _only_.


## Prerequisites

WSL2 requires that you are on Windows 11 or Windows 10 version 2004+ (Build 19041 and higher).
If you do not have at least these Windows versions, you must use [Lightsail for Research](../software-platforms/aws/index.md#lightsail-for-research-virtual-computing-with-aws) instead.


## Install WSL2

These instructions will provide you with both WSL2 and the Ubuntu app:

- WSL2 allows your computer to have a separate Linux subsystem
- Ubuntu is the default Linux distribution that gets installed when you install WSL2
    - The Ubuntu terminal application also acts as the Linux [terminal interface](../software-platforms/general-tools/using-the-terminal.md).
    - The Ubuntu terminal is the main way you will interact with the Linux side of your computer.
        - To learn more about the relationship between your Linux and Windows file systems, please refer to this article on [Working across Windows and Linux file systems](https://learn.microsoft.com/en-us/windows/wsl/filesystems).
        However, please be aware that for the purposes of OpenScPCA contribution, you should _only_ install files and perform analyses on the Linux file system.

To install WSL2, take the following steps:

1. In the Windows menu, search for the "Windows PowerShell" application.
Open it by clicking "Run as administrator".

1. Enter this command in PowerShell and hit enter:

    ```sh
    wsl --install
    ```

2. WSL2 will now install.
Along the way, you may get prompts asking if you allow the app to make changes to your device.
Always click "Yes" when you see these prompts.

1. Once WSL2 has finished installing, open the new Ubuntu app if it does not open automatically.
    - Ubuntu should prompt you to create a username and password which represent your WSL2 credentials.
      These are _independent_ of the username and password you already have set up on the Windows side of your computer.
      Changing one will not affect the other, but you can use the same username for both if you would like.
    - Note that when you type your password, no symbols will appear - this is expected!

    ??? question "Did Ubuntu not prompt you for a username and password?"
        If that prompt does not appear in Ubuntu, instead open PowerShell and run the command:

        ```sh
        wsl --unregister Ubuntu
        ```

        This will force Ubuntu to relaunch and prompt you for a username and password.

### Enable copy and paste

The standard `Ctrl+C/Ctrl+V` shortcuts for copy/paste do not natively work in Ubuntu.
To turn on copy/paste shortcuts, right-click on the Ubuntu window's title bar and select "Properties."

On the next screen, make sure `Use Ctrl+Shift+C/V as Copy/Paste` is checked on.
Now, you can use `Ctrl+Shift+C` for copy, and `Ctrl+Shift+V` for paste.
