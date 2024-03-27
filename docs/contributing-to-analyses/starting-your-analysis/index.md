# Starting Your Analysis

This section contains information about writing analysis modules.

For more information on using Git and GitHub while writing your analysis module, please see the [Working with Git](../working-with-git/index.md) documentation station.

### Determining and Managing Requirements for Your Analysis


When performing an analysis, there are two sets of requirements to be aware of and plan for:

1. **Software requirements**, often called dependencies, refer to the collection of software languages, packages, or libraries required for your analysis to run.
Tracking and reporting dependencies is important so that others can run your code and reproduce your results.

    - See the [Software Requirements](determining-software-requirements.md) page for approaches to track and manage dependencies.

1. **Compute requirements** refer to the computational resources, such as the number of CPUs or the amount of RAM, that are required for your analysis to run.
Gathering information about compute requirements helps us identify the computational resources we should use on Amazon Web Services.

    - See the [Compute Requirements](determining-compute-requirements.md) page for approaches to determine required computational resources.


### Creating your analysis module

We have provided a script to help you create your analysis modules with the expected file structure.
See the [Creating a module](./creating-a-module.md) page to learn how to use this script to make your analysis module folder.

### Modifying an analysis module

See the [Modifying an analysis module](STUB_LINK) page for specific considerations about steps you'll need to take to work with and modify the code in an existing module.

### Documenting an analysis module

A very important component of writing a module is to document that module.

Everyone benefits from documentation:

- It helps _you_ when you return to work on your module in the future
- It helps _other contributors_ get their bearings if they contribute to your module later
- It helps _the scientific community_ understand the scope, goals, and outputs of a module

See the [documenting your module](STUB_LINK) page to learn to document your OpenScPCA analysis module.
