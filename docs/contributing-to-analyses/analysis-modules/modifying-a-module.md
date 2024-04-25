# Modifying an existing module

When modifying or adding to an existing module, it's important to document all changes you make.

!!! note
    When you file an issue to track your work modifying an analysis module, you will use the [`Update an analysis` issue template](../../communications-tools/github-issues/issue-templates.md#update-an-analysis).

Please make sure to make the following updates when you work with an existing module:

- Fully document your changes [in the module's `README.md` file](./documenting-analysis.md).
This includes all of:
    - Changes to the behavior of any existing script or notebooks you modify
    - New scripts or notebooks you add to the module
    - Changes to which input data the module consumes
    - Changes to any output files (e.g., results or plots) the module generates
    - Changes to the module's [compute requirements](../determining-requirements/determining-compute-requirements.md)
    - Changes to instructions for setting up the module's run environment and/or running the module
- Add any new dependencies have been added to the module's [software environment](../determining-requirements/determining-software-requirements.md)
- Update current documentation within any scripts or notebooks that you modify so that it describes code's updated behavior.
