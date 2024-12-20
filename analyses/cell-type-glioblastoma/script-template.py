#!/usr/bin/env python3

# cell-type-glioblastoma Analysis Script
# Author Name
# Date
#

# This file is a template for an analysis script written in Python.
#
# It is a good idea to start with a brief introduction to the contents of the
# script, including the purpose of the analysis, the data used, and the methods
# applied.
#
# Please also include examples of of how the script should be invoked, with any
# required arguments or options.
#
# Replace this set of comment with your own introduction, and be sure to update
# the Title, Author Name, and Date at the top of the document.
# Don't forget to rename this file as well!
#
# Include comments in your code to explain your steps and to document any
# assumptions you are making.

# Load modules
#
# Load required Python modules at the top of your script
# We have included the standard `pathlib` module and the `session_info` module
# that we will be using at the bottom of this notebook to record the versions of
# the modules used in this analysis.
#
# Do not install modules here; only load them with `import` statements.
# Avoid renaming modules with `as` statements, unless you are performing a
# standard renaming (e.g., `import pandas as pd`).

import contextlib
import pathlib

import session_info

# Define functions


def find_git_root():
    """Find the root directory of the git repository from the calling location of this script"""
    repo_root = pathlib.Path.cwd()
    while not (repo_root / ".git").is_dir():  # search for the .git directory
        repo_root = repo_root.parent
        if repo_root == "/":
            raise FileNotFoundError("Could not find the repository root directory")
    return repo_root


# Below we set some standard paths that may be useful for your script
# find the repository root directory
repo_root = find_git_root()

# set module path (using pathlib)
module_root = repo_root / "analyses" / "cell-type-glioblastoma"

# set current data directory
data_dir = repo_root / "data" / "current"

# set results and plots directories (using the analysis project file to find root)
results_dir = module_root / "results"
plots_dir = module_root / "plots"

# Set other paths as needed:

# Input files

# Output files
session_info_path = results_dir / "session_info.txt"

# Main body of the script goes here


# As the last step, record the versions of the modules and dependencies
# that were used in this analysis
with open(session_info_path, "w") as f:
    with contextlib.redirect_stdout(f):  # direct the session_info output to a file
        session_info.show(dependencies=True)
