# Release notes for OpenScPCA-analysis

This document contains release notes for versioned releases of the OpenScPCA-analysis repository.
We will also periodically release and tag dated snapshots of the repository, but these will not be reflected in this document.

<!--
Add new release notes in reverse numerical order (newest first) below this comment

You may want to add temporary notes here for tracking as features are added, before a new release is ready.
-->


## v0.1.0

This is the initial versioned release of the OpenScPCA-analysis repository.
The repository at this stage should be generally complete with respect to infrastructure, but the analysis is still only in the early stages.

With respect to infrastructure, the repository contains the following components:

- detailed documentation in the `docs` directory explaining how to interact with the OpenScPCA project and how to set up and run analyses
- a `create-analysis-module.py` script for setting up new analysis modules
- `download-data.py` and `download-results.py` scripts to download data and results from the OpenScPCA project
- template notebooks, scripts, environment files, and Docker images for analysis modules

The repository contains the following analysis modules that are now in a relatively mature state:

- hello-R
- hello-python
- simulate-sce
- doublet-detection

The following additional modules are in active development:

- cell-type-ewings
- metacells
