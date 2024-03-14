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

The instructions provided here are suitable for the vast majority of contributors, but please refer to [this section](#advanced-considerations) for advanced considerations about R installation.


### Install R on macOS

1. Navigate to the [macOS download page on the CRAN website](https://cran.r-project.org/bin/macosx/).

2. Download the R package that matches your computer's architecture, and follow all installation instructions.
    - If you're on an Apple silicon (M1-3) Mac, install R from the link `R-X.Y.Z-arm64.pkg`, where `X.Y.Z` is the specific R version.
    - If you're on an Intel Mac, install R from the link `R-X.Y.Z-x86_64.pkg`, where `X.Y.Z` is the specific R version.


Next, you need to install additional tools that will let you build R packages that you are likely to encounter in OpenScPCA.
These tools are available from the [R for macOS website](https://mac.r-project.org/tools/) and are listed under **Mandatory Tools**.

1. Install "**Xcode** developer tools from Apple" as follows:
    - Launch a Terminal window by searching for "terminal" in spotlight
      <figure markdown="span">
        ![Launch the Terminal application.](../../img/terminal-spotlight.png){width="425"}
      </figure>
    - Enter the command `sudo xcode-select --install` into Terminal and hit enter
    - You will then be prompted to enter your computer's password
  When you type, no symbols will appear - this is expected!
  Hit enter after your password, and the developer tools will proceed to install
        - If you receive a message that begins with `xcode-select: note: Command line tools are already installed.`, then you have previously installed this on your system and are all set.

1. Install the **GNU Fortran compiler**.
Click the link to download the installer package and follow all installation instructions.


### Install R on Windows


1. Navigate to the [Windows download page on the CRAN website](https://cran.r-project.org/bin/windows/base/), and follow instructions to download and install R.

2. Next, you need to install Rtools, which provides additional tools that will let you build R packages that you are likely to encounter in OpenScPCA.
    - Navigate to the [Rtools download page](https://cran.r-project.org/bin/windows/Rtools/)
    - Click the Rtools version that matches the R version you just downloaded
    - On the next page, click the link `RtoolsXY installer` (where `XY` is the specific version you clicked) to download Rtools, and follow all installation instructions




## Install the RStudio IDE

If you do not already have RStudio installed, you need to [install the free RStudio Desktop](https://posit.co/download/rstudio-desktop/).
This website shows Step 1 as installing R itself, which you can skip and proceed straight to Step 2, Install RStudio.

Click the download link and follow all instructions to complete the installation.

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

If you need to have multiple versions of installed on your computer, we recommend using the [`rig`](https://github.com/r-lib/rig) R Installation Manager.