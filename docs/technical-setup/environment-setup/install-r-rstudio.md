# Installations for development with R

This page provides instructions on what you need to install to be able to contribute analysis modules that use R.

You will need to have the following items installed:

- R version `4.3.3`
- The RStudio IDE
- The `renv` R package

Please read these instructions carefully even if you already have any of these items installed on your system because there are a few additional items you need to make sure you have.


## Install R

`OpenScPCA` will use R version `4.3.3`.
This means that during your time contributing to `OpenScPCA`, you should not plan to update your R installation with a future release (however, please see [this section](#advanced-considerations) if you need to maintain multiple versions of R on your system).



### Install R on macOS

1. From the [CRAN website](https://cran.r-project.org/), click the link "Download R for macOS".

1. Scroll to the page section titled `R 4.3.3 "Angel Food Cake" released on 2024/02/29`. **TODO: What happens when this is no longer the latest release?**

1. Download the R package that matches your computer, and follow all installation instructions.
    - If you're on an Apple silicon (M1-3) Mac, install R from `R-4.3.3-arm64.pkg`
    - If you're on an Intel Mac, install R from `R-4.3.3-x86_64.pkg`

1. Next, you will need to install additional tools that will let you build R packages that you are likely to encounter in `OpenScPCA`.
These tools are available from the [Mac R Project](https://mac.r-project.org/tools/) and are listed under **Mandatory Tools**.
    - First, install "**Xcode** developer tools from Apple" as follows:
      - Launch a `Terminal` window by searching for "terminal" in spotlight
        <figure markdown="span">
            ![Launch the Terminal application.](../../img/terminal-spotlight.png){width="425"}
        </figure>
      - Enter the command `xcode-select --install` into Terminal and hit enter.
    - Next, install the **GNU Fortran compiler**.
    Click the link to download the installer package and follow all installation instructions.


### Install R on Windows

## Advanced considerations

If you need to have multiple versions of installed on your computer, we recommend using the [`rig`](https://github.com/r-lib/rig) R Installation Manager.