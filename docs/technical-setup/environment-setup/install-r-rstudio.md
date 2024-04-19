# Install R, RStudio, and other dependencies

This page provides instructions on the tools you need to install to be able to contribute to or explore analysis modules with R.
If you do not plan on using R as part of your OpenScPCA contribution, you do not need these tools.

If you plan to use R, you will need to have the following tools installed:

- R
- The RStudio IDE
- The `renv` and `BiocManager` R packages

Please read these instructions carefully even if you already have any of these installed on your system because there are a few additional tools you need to make sure you have.

## Install R

The official OpenScPCA project version of R is `4.3.3`, so we recommend that you have at least this version installed.
Please be aware that if you install a different R version and run existing modules, you may get different results.

The instructions provided here are suitable for the vast majority of contributors, but please refer to our section on [advanced considerations about R installation](#advanced-considerations) for more information.


### Install R on macOS

1. Navigate to the [macOS download page on the CRAN website](https://cran.r-project.org/bin/macosx/).

2. Download the R package that matches your computer's architecture, and follow all installation instructions.
    - If you're on an Apple silicon (M1-3) Mac, install R from the link `R-X.Y.Z-arm64.pkg`, where `X.Y.Z` is the specific R version.
    - If you're on an Intel Mac, install R from the link `R-X.Y.Z-x86_64.pkg`, where `X.Y.Z` is the specific R version.


Next, you need to install additional tools that will let you build R packages that you are likely to encounter in OpenScPCA.
These tools are available from the [R for macOS website](https://mac.r-project.org/tools/) and are listed under **Mandatory Tools**.

1. Install "**Xcode** developer tools from Apple" as follows:
    - Launch a [Terminal window](../../software-platforms/general-tools/using-the-terminal.md)
    - Enter the command `sudo xcode-select --install` into Terminal and hit enter
    - You will then be prompted to enter your computer's password
  When you type, no symbols will appear - this is expected!
  Hit enter after your password, and the developer tools will proceed to install
        - If you receive a message that begins with `xcode-select: note: Command line tools are already installed.`, then you have previously installed this on your system and are all set with this step.

1. Install the **GNU Fortran compiler**.
Click the link on the [Tools - R for Mac OS](https://mac.r-project.org/tools/) page to download the `gfortran` installer package and follow all installation instructions.

### Install R on Windows with WSL2

These instructions will install R on the WSL2 side of your computer.

1. Copy and paste this code into the Ubuntu terminal to install R:

    ```sh
    # install a few helper packages needed to install R and certain R packages
    sudo apt install --no-install-recommends software-properties-common dirmngr libssl-dev

    # add the signing key (by Michael Rutter) for these repos
    wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | sudo tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc

    # add the R 4.0 repo from CRAN
    sudo add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"

    # install R
    sudo apt install --no-install-recommends r-base r-base-dev
    ```

1. To confirm R was successfully installed, run `R` in the Ubuntu terminal.
This should launch the R console in terminal; hit `Ctrl+D` to quit.
    - If you instead get an error message that the command was not found, your R installation was not successful; please feel free to use the [Troubleshooting GitHub Discussions](../../communications-tools/index.md#ask-questions) to get some help!

1. Finally, copy and paste this line into the Ubuntu terminal.
    ```sh
    echo "options(repos = list(CRAN='https://p3m.dev/cran/__linux__/$(lsb_release -cs)/latest'))" >> ~/.Rprofile
    ```
      - This line sets the default R package repository to P3M, the Posit Public Package Manager, instead of the default of CRAN in your [`.Rprofile`](https://support.posit.co/hc/en-us/articles/360047157094-Managing-R-with-Rprofile-Renviron-Rprofile-site-Renviron-site-rsession-conf-and-repos-conf) file.
       - It will _dramatically_ streamline R package installation by providing you with pre-built package binaries, removing the need to install lots of additional system library dependencies on your computer.



## Install the RStudio IDE

### Install the RStudio IDE on macOS

If you do not already have RStudio installed, you need to [install the free RStudio Desktop](https://posit.co/download/rstudio-desktop/).
This website shows Step 1 as installing R itself, which you can skip and proceed straight to Step 2, Install RStudio.

Click the download link and follow all instructions to complete the installation.


### Install RStudio Server on Windows with WSL2

Because there is no native way to use the RStudio IDE within WSL2's Ubuntu operating system, you will instead need to install the RStudio Server.

This will provide essentially the same experience as working with RStudio Desktop, except you will access the IDE through your browser instead of through a built-in GUI.

To install RStudio Server, run these commands in the Ubuntu terminal:

1. Download the RStudio Server installer file:

    ```sh
    wget "https://rstudio.org/download/latest/stable/server/$(lsb_release -cs)/rstudio-server-latest-amd64.deb"
    ```

1. Install RStudio Server:

    ```
    sudo apt install ./rstudio-server-latest-amd64.deb
    ```

    Note that this line will also start the server automatically when the installation is complete.

1. You can now safely remove the installer file:


    ```sh
    rm rstudio-server-latest-amd64.deb
    ```


#### Using the RStudio Server

To access the RStudio Server, the server must be running.
You can start and stop the server by running the following lines in the Ubuntu terminal.

- Start RStudio Server

    ```sh
    sudo rstudio-server start
    ```

- Stop RStudio Server

    ```sh
    sudo rstudio-server stop
    ```

To use the RStudio Server, navigate to `http://localhost:8787` in your browser.
You will be prompted to sign into the server with a username and password; use the username and password credentials you created when you [installed Ubuntu](../install-wsl2.md#installation-instructions).

Once you log in, you will see an RStudio IDE in the browser which you can use as you normally would use RStudio Desktop!
Because the server is installed on the WSL2 side of your computer, you can fully interact with your Linux file system via RStudio Server's built-in terminal just the same as the Ubuntu terminal application.




## Install R package dependencies

There are two R packages that OpenScPCA frequently uses that we recommend installing before getting started:

- `renv`
    - OpenScPCA uses [the `renv` package](https://rstudio.github.io/renv/articles/renv.html) to manage R packages in each analysis module.
- `BiocManager`
    - You need to use [Bioconductor packages](https://bioconductor.org/) to interact with ScPCA data in R environments.
  [The `BiocManager` package](https://cran.r-project.org/web/packages/BiocManager/vignettes/BiocManager.html) is used to install Bioconductor packages.


To install these R packages, open RStudio and enter the following command into the Console:

```
install.packages(c("renv", "BiocManager"))
```


## Advanced considerations

If you need to have multiple versions of R installed on your computer, we recommend using the [`rig` R Installation Manager](https://github.com/r-lib/rig).
