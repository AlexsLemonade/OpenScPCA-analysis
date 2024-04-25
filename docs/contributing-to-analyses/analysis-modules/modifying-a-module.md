# Modifying an existing module


If you are modifying a module that already exists, then you will first need to file an issue using the [`Update an analysis` issue template](../../communications-tools/github-issues/issue-templates.md#update-an-analysis). 
Then you can create a new branch and commit any code changes. 

Be sure that any new code changes are accompanied by updates to documentation. 
You will also want to capture any new software dependencies. 

The following changes should be included in any PR's that modify existing modules: 

- Fully document your changes [in the module's `README.md` file](./documenting-analysis.md).
This includes all of:
    - Changes to the behavior of any existing script or notebooks you modify
    - New scripts or notebooks you add to the module
    - Changes to which input data the module consumes
    - Changes to any output files (e.g., results or plots) the module generates
    - Changes to the module's [compute requirements](../determining-requirements/determining-compute-requirements.md)
    - Changes to instructions for setting up the module's run environment and/or running the module
- Make sure any new dependencies have been added to the module's [software environment](../determining-requirements/determining-software-requirements.md)
- Either update or add `README.md` files within any sub-folders present in the module. 
For example, if you added a new script to the `scripts` directory, be sure to add a `scripts/README.md` that describes the usage of the new script. 
