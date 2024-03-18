# Determining Compute Requirements

We use the term **compute requirements** to refer to the computational resources, such as the number of CPUs or the amount of RAM, that are required for your analysis to run.
Reporting the compute requirements for running your analysis to project organizers helps us plan for allocating resources on Amazon Web Services (AWS) and to determine how to review [your pull requests](STUB_LINK PRs docs) (e.g., using our laptops vs. using AWS).

We will prompt you to report compute requirements throughout the project as described [below](#reporting-compute-requirements).

## Determining compute requirements

### Prior success or failure

A great place to start with determining compute requirements is what computer(s) you can successfully run your analysis on.

* Have you run this kind of analysis successfully on your laptop?
  * You can use your laptop's specs when reporting compute requirements.
* Have you successfully run a similar analysis before on your institution's High Performance Computing (HPC) cluster?
  * Report the characteristics of the kind of node you typically use.

Even failure can be very helpful!
For instance, maybe you ran out of memory while running this analysis on a particular computer.
The amount of RAM on that computer is helpful information, because it tells us that the analysis probably requires more RAM.

### Read the docs

If you know you will be using a particular software package for your analysis, review the package's documentation to see if the authors and maintainers give you any guidance about the compute requirements for running that package.

For example, the [`scvi-tools` documentation](https://docs.scvi-tools.org/en/stable/index.html) tells you that [using a GPU helps the analysis run faster](https://docs.scvi-tools.org/en/stable/installation.html#gpu).
Methods or benchmarking papers in the literature may also be a good resource for figuring out compute requirements.

### How to monitor requirements

#### `htop`

To get an overview of computational resource usage when running an analysis, you can use [`htop`](https://htop.dev/) on Linux and macOS.

`htop` will be already installed on machines run via [Lightsail for Research](STUB_LINK), but you can [install it via conda](https://anaconda.org/conda-forge/htop) (see our [documentation on conda](software-requirements.md#adding-software-to-the-environment-and-tracking-installed-software)).

Check out [this beginner's guide to `htop`](https://spin.atomicobject.com/htop-guide/) for an introduction.

### Mind the file sizes

Although we are most concerned with finding the right size computer to run an analysis, we also need to be concerned about data storage, particularly when it comes to [volumes on AWS](STUB_LINK volumes on AWS).

We generally will not ask you to report input file sizes because we can calculate them if we know what project(s) or sample(s) you are analyzing!

However, if your module code generates large intermediate files, document that information in the README of your module to help others who want to run the code.

Check out [this blog post on finding file sizes in Linux](https://monovm.com/blog/how-to-find-the-file-size-in-linux/) for a variety of tips on calculating file and directory sizes.

## Reporting compute requirements

### Discussions, issues, and pull requests

We use templates throughout the project to prompt project maintainers and contributors to provide helpful information.
Here, we describe when and where to document compute requirements on GitHub discussions, issues, and pull requests.

Please be sure to provide this information when prompted.

* [Proposing a new analysis in Discussions?](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/new?category=propose-a-new-analysis)
  * Report what you know about computational resources required in the "Other Details" field
* [Proposing to modify an existing analysis in Discussions?](https://github.com/AlexsLemonade/OpenScPCA-analysis/discussions/new?category=modify-an-existing-analysis)
  * Report what you know about computational resources required in the "Other Details" field
* [Filing an issue to start an analysis?](STUB_LINK)
  * Report what you know about computational resources required in the "What computational resources will you require?" field
* [Filing an issue to modify an analysis?](STUB_LINK)
  * Report what you know about computational resources required in the "Will you require different computational resources beyond what the analysis module already uses?" field
* [Filing a pull request?](STUB_LINK)
  * Fill out the section explaining to reviewers what computational resources are required to run the code ("What are the computational requirements to be able to run the code in this PR?").

You can always add a comment to your existing discussion post or issue about compute requirements if you have gathered more information since filing!

### READMEs

Every analysis module in the project includes its own README file.
Document the compute requirements in the README of your module. (STUB_LINK)
