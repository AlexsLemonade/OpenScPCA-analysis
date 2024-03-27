# Documenting your analysis

Each analysis module should include a `README.md` file with important information about the analysis and instructions for running the analysis module.
This file can be particularly helpful to other contributors who want to use your analysis module and to you when you revisit a module you've written in the past.

We have provided a [template `README.md`](../../../templates/README.md) that will be automatically included when you [create a new analysis module](./creating-a-module.md).
You should fill out the template for your analysis module and provide the following information for each analysis:

- Description
- Usage (How to run the analysis)
- Input files
- Output files
- Software requirements
- Computational resources

We recommend that as you continue to work on the analysis, you update the `README.md` file with each change.
For example, if you add a new step in the analysis module, your [pull request](../creating-pull-requests/index.md) should include the code changes and an update to the `README.md` containing a description of that step.

!!! note "Is this your first analysis module?"
    If you are creating your first analysis module, we recommend that you [scope your first pull request](../creating-pull-requests/scoping-pull-requests.md) to contain only:

    - The analysis module skeleton [created by running `create-analysis-module.py`](./creating-a-module.md)
    - A `README.md` file containing at minimum, a description of the analysis
        - Other headers in the `README.md`, like `Usage`, can be filled out as the analysis proceeds
