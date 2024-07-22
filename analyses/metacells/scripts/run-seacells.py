#!/usr/bin/env python3

"""
Seacells Analysis Script
Joshua Shapiro
2024-07-18

This script runs the SEACell algorithm on the given dataset.
"""

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

import argparse
import contextlib
import pathlib

import anndata
import SEACells
import session_info


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the SEACell algorithm on the given dataset."
    )
    parser.add_argument(
        "datafile", type=pathlib.Path, help="The input data in H5AD format."
    )
    parser.add_argument("output", type=pathlib.Path, help="The output file path.")
    parser.add_argument("logfile", type=pathlib.Path, help="File path for log outputs")

    args = parser.parse_args()

    adata = AnnData.read(args.datafile)

    # As the last step, record the versions of the modules and dependencies
    # that were used in this analysis
    with open(args.logfile, "w") as f:
        with contextlib.redirect_stdout(f):  # direct the session_info output to a file
            session_info.show(dependencies=True)


if __name__ == "__main__":
    main()
